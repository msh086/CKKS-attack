package main

import (
	"flag"
	"fmt"
	"github.com/ldsec/lattigo/v2/ckks"
	"math"
	"math/bits"
	"math/cmplx"
	"math/rand"
	"time"
)

var modeFlag string
var lognFlag uint64
var scaleFlag float64
var logScaleFlag uint64
var funcTypeFlag string

var maclaurinSeries = map[string][]complex128{
	"exponent": {
		complex(1.0, 0),
		complex(1.0, 0),
		complex(1.0/2, 0),
		complex(1.0/6, 0),
		complex(1.0/24, 0),
		complex(1.0/120, 0),
		complex(1.0/720, 0),
		complex(1.0/5040, 0),
		complex(1.0/40320, 0)},
		//1./2, 1./4, 0, -1./48, 0, 1./480, 0, -17./80640, 0, 31./1451520, 0
	"sigmoid": {
		complex(1.0/2, 0),
		complex(1.0/4, 0),
		complex(0, 0),
		complex(-1.0/48, 0),
		complex(0, 0),
		complex(1.0/480, 0),
		complex(0, 0),
		complex(-17.0/80640, 0),
		complex(0, 0),
		complex(31.0/1451520, 0),
	},
}

func init() {
	flag.StringVar(&modeFlag, "mode", "attack", "mode to run, available options: attack|defense")
	flag.Uint64Var(&lognFlag, "logn", 16, "log2(N), N is the ring dimension")
	flag.Float64Var(&scaleFlag, "scale", 1 << 40, "scale for encoding")
	flag.Uint64Var(&logScaleFlag, "logscale", 0, "log2(scale), has higher priority than '-scale', " +
		"when set to 0, the value of '-scale' will be used instead")
	flag.StringVar(&funcTypeFlag, "func", "", "function to be evaluated. {exponent, sigmoid, variance}")
}

func printFLags() {
	fmt.Printf("mode = %s, logn = %d, scale = %d, log scale = %d", modeFlag, lognFlag, scaleFlag, logScaleFlag)
}

func checkFlags() {
	if modeFlag != "attack" && modeFlag != "defense" {
		panic(fmt.Sprintf("unknown mode %s", modeFlag))
	}
	if lognFlag <= 1 {
		panic(fmt.Sprintf("%d is too small for logn", lognFlag))
	}
	if lognFlag > 16 {
		panic(fmt.Sprintf("%d is too large for logn", lognFlag))
	}
	if logScaleFlag != 0 {
		scaleFlag = float64(uint64(1) << logScaleFlag)
	}
}

func infNorm(arr []complex128) float64 {
	var max float64 = 0
	for _, ele := range arr {
		im, re := imag(ele), real(ele)
		sqrL2 := im * im + re * re
		if max < sqrL2 {
			max = sqrL2
		}
	}
	return math.Sqrt(max)
}

func homoMean(ciphertext *ckks.Ciphertext, evaluator ckks.Evaluator, rotk *ckks.RotationKeys) *ckks.Ciphertext {
	ctxtCopy := ciphertext.CopyNew().Ciphertext()
	tmp := ciphertext.CopyNew().Ciphertext()
	logn := bits.Len(uint(len(ctxtCopy.Value()[0].Coeffs[0]))) - 1
	for step := uint64(1); logn > 0; logn, step = logn - 1, step << 1  {
		evaluator.Rotate(ctxtCopy, step, rotk, tmp)
		evaluator.Add(ctxtCopy, tmp, ctxtCopy)
	}
	// NOTE: passing 1.0/len(...) will get an int of value 0 !!!
	evaluator.MultByConst(ctxtCopy, 1.0/float64(len(ciphertext.Value()[0].Coeffs[0])), ctxtCopy)
	evaluator.RescaleMany(ctxtCopy, 1, ctxtCopy)
	return ctxtCopy
}

func homoVar(ciphertext *ckks.Ciphertext, evaluator ckks.Evaluator,
		evk *ckks.EvaluationKey, rotk *ckks.RotationKeys) *ckks.Ciphertext {
	ctxtCopy := evaluator.ConjugateNew(ciphertext, rotk)
	evaluator.MulRelin(ciphertext, ctxtCopy, evk, ctxtCopy)
	if err := evaluator.RescaleMany(ctxtCopy, 1, ctxtCopy); err != nil {
		panic(err)
	}
	meanSquare := homoMean(ctxtCopy, evaluator, rotk) // E(|x|^2)
	mean := homoMean(ciphertext, evaluator, rotk)
	evaluator.Conjugate(mean, rotk, ctxtCopy)
	evaluator.MulRelin(mean, ctxtCopy, evk, mean) // |E(x)|^2
	if err := evaluator.RescaleMany(mean, 1, mean); err != nil {
		panic(err)
	}
	evaluator.Sub(meanSquare, mean, meanSquare)
	return meanSquare
}

func sampleUnitCircle() complex128 {
	cdfSamp := float64(rand.Uint64()) / math.MaxUint64
	rSamp := math.Sqrt(cdfSamp)
	rhoSamp := float64(rand.Uint64()) / math.MaxUint64 * 2 * math.Pi
	return cmplx.Rect(rSamp, rhoSamp)
}

func evalPointPoly(val complex128, coeffs []complex128) (res complex128) {
	res = coeffs[0]
	pow := val
	for _, coeff := range coeffs[1:] {
		res += pow * coeff
		pow *= val
	}
	return
}

func evalVecPoly(vec []complex128, coeffs []complex128) {
	for idx := range vec {
		vec[idx] = evalPointPoly(vec[idx], coeffs)
	}
}

func evalVecVar(vec []complex128) (Var float64) {
	var mean complex128 = 0
	Var = 0
	for _, val := range vec {
		mean += val
		Var += imag(val) * imag(val) + real(val) * real(val)
	}
	mean /= complex(float64(len(vec)), 0) // E(x)
	Var /= float64(len(vec)) // E(|x|^2)
	Var -= imag(mean) * imag(mean) + real(mean) * real(mean) // E(|x|^2) - |E(x)|^2
	return
}

//func evalPolySimple(src *ckks.Ciphertext, dst *ckks.Ciphertext, coeffs []complex128,
//		evaluator ckks.Evaluator, encoder ckks.Encoder, relinKey *ckks.EvaluationKey, parameters *ckks.Parameters) {
//	polyDeg := len(coeffs)
//	bitsLen := bits.Len64(uint64(polyDeg))
//	powerOf2s := make([]*ckks.Ciphertext, bitsLen)
//	powerOf2s[0] = src
//	for i := 1; i < bitsLen; i++ {
//		powerOf2s[i] = evaluator.MulRelinNew(powerOf2s[i - 1], powerOf2s[i - 1], relinKey)
//		evaluator.Rescale()
//	}
//	dst = ckks.NewCiphertext(parameters, src.Degree(), src.Level(), src.Scale())
//	dst.Copy(src.El())
//	tmp := ckks.NewPlaintext(parameters, src.Level(), src.Scale())
//	valueVec := make([]complex128, 1 << parameters.LogSlots())
//	setVec := func (val complex128) {
//		for idx := range valueVec {
//			valueVec[idx] = val
//		}
//	}
//	setVec(coeffs[1])
//	encoder.Encode(tmp, valueVec, parameters.LogSlots())
//	evaluator.MulRelin(dst, tmp, relinKey)
//}

func example() {
	flag.Parse()

	printFLags()

	checkFlags()

	var start time.Time
	var err error

	LogN := lognFlag //uint64(16)
	LogSlots := lognFlag - 1 //uint64(15)

	LogModuli := ckks.LogModuli{
		LogQi: []uint64{55, 40, 40, 40, 40, 40, 40, 40},
		LogPi: []uint64{45, 45},
	}

	Scale := scaleFlag //float64(1 << 40)

	// Note:
	//  although NTT only requires that q mod n = 1, but here primes are chosen s.t. q mod 2n = 1
	//  shijie provides a possible explanation for this:
	//	 the ring R is Z[x]/(x^n+1), but we expand to Z[x]/(x^2n-1) and calculate NTT there
	//	 since x^2n-1 is the typical form for 'primitive roots of unity', while x^n+1 is not
	//
	//  ref to ring.go line 117 genNTTParams
	params, err := ckks.NewParametersFromLogModuli(LogN, &LogModuli)
	if err != nil {
		panic(err)
	}
	params.SetScale(Scale)
	params.SetLogSlots(LogSlots)

	fmt.Println()
	fmt.Println("=========================================")
	fmt.Println("         INSTANTIATING SCHEME            ")
	fmt.Println("=========================================")
	fmt.Println()

	start = time.Now()

	/*
		Why I am looking into all kinds of modular reduction:
			the MForm behaviour differs from what I read on websites, where I'm told both multiplication operands
			need to be transformed into Montgomery form to perform multiplication, then a reduction and a real
			mod operation take place, getting the final result
			But in lattigo only either operand is turned into Montgomery form
			below is a detailed explanation of modular reductions in lattigo

		MRedParams: q^-1 mod 2^64 (Note: the k = (r(r^-1 mod q) - 1) / q in Montgomery is q^-1 mod r)
		BRedParams: floor(2^128/q) stored in 2 uint64s
		MForm: using Barret Reduction to compute (a * 2^64 mod q), here:
			Barret k = 64, modulus = q
			Montgomery: r = 2^64, modulus = q
		MRed: takes x, y, q, qinv, output r
			x is the output of MForm(a, q), i.e. x = (a * 2^64 mod q), which is already in Montgomery form
			y is NOT in Montgomery form
			first compute double-precision integer m := x * y \in [0, q^2)
			then compute double-precision h := (x * qinv mod 2^64) * q
			according to https://www.nayuki.io/page/montgomery-reduction-algorithm, m and h are summed up,
			then divided by 2^64 to get r
			but here 2*2^64*q may exceed 2^128(really? see Note below) and it involves double-precision integer addition
			so instead q + (m - h) / 2^64 is computed, with the result in range (0, 2q)
			then a simple subtraction is performed
			// Note: I think q should be less than 2^63, or reduction from (0, 2q) into q won't work...
			// 	so actually 2*2^64*q won't exceed 2^128
	*/
	kgen := ckks.NewKeyGenerator(params)

	sk := kgen.GenSecretKey()

	rlk := kgen.GenRelinKey(sk)

	rotk := kgen.GenRotationKeysPow2(sk)
	kgen.GenRotationKey(ckks.Conjugate, sk, 1, rotk)

	encryptor := ckks.NewEncryptorFromSk(params, sk)

	decryptor := ckks.NewDecryptor(params, sk)

	encoder := ckks.NewEncoder(params)

	evaluator := ckks.NewEvaluator(params)

	fmt.Printf("Done in %s \n", time.Since(start))

	fmt.Println()
	fmt.Printf("CKKS parameters: logN = %d, logSlots = %d, logQP = %d, levels = %d, scale= %f, sigma = %f \n", params.LogN(), params.LogSlots(), params.LogQP(), params.MaxLevel()+1, params.Scale(), params.Sigma())

	fmt.Println()
	fmt.Println("=========================================")
	fmt.Println("           PLAINTEXT CREATION            ")
	fmt.Println("=========================================")
	fmt.Println()

	start = time.Now()

	//r := float64(16)
	//pi := 3.141592653589793

	slots := params.Slots()

	values := make([]complex128, slots)
	for i := range values {
		//values[i] = complex(2*pi, 0)
		values[i] = sampleUnitCircle()
	}

	/*
	   a plaintext polynomial is in NTT form, it's represented in a (nbModuli x N) matrix of uint64
	*/
	plaintext := ckks.NewPlaintext(params, params.MaxLevel(), params.Scale())
	/*
		//in call sequence
		//	ckks/encoder.go line 132 encoder::Encode
		//	ckks/encoder.go line 186 encoder::ScaleUp
		//	ckks/utis.go line 38 scaleUpVecExact
		//		if the scaled value n*values[i] exceed pow(2.0, 64), multi-precision integers are used
		//		to prevent overflow during conversion from float64 to uint64
		//		the logic should be "if abs(n * values[i]) > pow(2.0, 64)",
		//		but the code is actually "if n * values[i] > pow(2.0, 64)",
		//		could this possibly be a bug?
		NOTE: I sent an issue #102 concerning this bug, now it's fixed in the coming commit
	*/
	encoder.Encode(plaintext, values, params.LogSlots())

	fmt.Printf("Done in %s \n", time.Since(start))

	fmt.Println()
	fmt.Println("=========================================")
	fmt.Println("              ENCRYPTION                 ")
	fmt.Println("=========================================")
	fmt.Println()

	start = time.Now()

	/*
		here the encryptor is an instance of skEncryptor, calling EncryptNew will generate a new random element
		of R_q in NTT form
		the modulus are deliberately chosen s.t. q_i = 1 mod 2N
		i.e. NTT modulus is the same as non-NTT modulus
		which means sampling a random element in NTT form is essentially the same as sampling a random element
		in non-NTT form
	*/
	ciphertext := encryptor.EncryptNew(plaintext)

	fmt.Printf("Done in %s \n", time.Since(start))

	printDebug(params, ciphertext, values, decryptor, encoder)

	fmt.Println()
	fmt.Println("===============================================")
	fmt.Printf("       EVALUATION of func on %d values\n", slots)
	fmt.Println("===============================================")
	fmt.Println()

	start = time.Now()

	coeffs, ok := maclaurinSeries[funcTypeFlag]
	if ok {
		poly := ckks.NewPoly(coeffs)

		// TODO: evaluator.EvaluatePoly is not so simple
		//  having something to do with 'Baby-step, Giant-step Algorithm'
		//  I know BSGS algorithm in computing discrete logarithm
		//  and I can vaguely guess what they want to do
		//  but the code is so complicated for me to understand
		//  possible ref: https://specfun.inria.fr/bostan/publications/exposeJNCF.pdf page 78
		//  Paterson-Stockmeyer algorithm
		//
		// FIXME: another little bug (not necessarily) there...
		//  in polynomial_evaluation.go line 270, recurse
		//  shift by (logSplit - 1) when logSplit == 1 will cause a panic
		if ciphertext, err = evaluator.EvaluatePoly(ciphertext, poly, rlk); err != nil {
			panic(err)
		}
	} else if funcTypeFlag == "variance" {
		ciphertext = homoVar(ciphertext, evaluator, rlk, rotk)
	}

	fmt.Printf("Done in %s \n", time.Since(start))

	if ok {
		evalVecPoly(values, coeffs)
	} else if funcTypeFlag == "variance" {
		plainVar := evalVecVar(values)
		for idx := range values {
			values[idx] = complex(plainVar, 0)
		}
	}

	printDebug(params, ciphertext, values, decryptor, encoder)

	fmt.Println()
	fmt.Println("=========================================")
	fmt.Println("         DECRYPTION & DECODING           ")
	fmt.Println("=========================================")
	fmt.Println()

	start = time.Now()

	outPlain := decryptor.DecryptNew(ciphertext)
	outValues := encoder.Decode(outPlain, params.LogSlots())

	fmt.Printf("Done in %s \n", time.Since(start))

	fmt.Println("this is the precision stats between expected & actual values")
	printDebug(params, ciphertext, values, decryptor, encoder)

	// now attack
	fmt.Println("log2(scale) = ", math.Log2(ciphertext.Scale()))
	fmt.Println("log2(inf norm) = ", infNorm(outValues))
	ptxt1 := ckks.NewPlaintext(params, ciphertext.Level(), ciphertext.Scale())
	ptxt2 := ckks.NewPlaintext(params, ciphertext.Level(), ciphertext.Scale())

	if modeFlag == "defense" {
		fmt.Println("adding extra error for defense...")

		expEnc := ckks.CastExposedEncryptor(encryptor)
		expEnc.AddGaussianNoise(outPlain.El())

		noisy := encoder.Decode(outPlain, params.LogSlots())
		fmt.Println("this is the precision stats between actual & noisy values")
		printDebugBrief(params, outValues, noisy)
		fmt.Println("this is the precision stats between expected & noisy values")
		printDebugBrief(params, values, noisy)
		outValues = noisy
	}

	encoder.EncodeNTT(ptxt1, outValues, params.LogSlots())
	fmt.Println("re-encoding ok? ", ckks.IsSame(ptxt1.El(), outPlain.El()))
	// now ptxt1 = m+e
	exposed := ckks.CastExposedEvaluator(evaluator)
	cipherB := ckks.GetCiphertextPartAsElement(ciphertext, 0)
	exposed.SubPlain(ptxt1.El(), cipherB, ptxt1.El())
	// now ptxt1 = m+e-b
	cipherA := ckks.GetCiphertextPartAsElement(ciphertext, 1)
	if err := exposed.ModInverse(cipherA, ptxt2.El()); err != nil {
		fmt.Println(err)
		return
	}
	// now ptxt2 = a^-1
	exposed.MultElement(ptxt2.El(), ptxt1.El(), ptxt1.El())
	// now ptxt1 = a^-1*(m+e-b)
	exposed.MFormElement(ptxt1.El(), ptxt1.El())
	skEle := ckks.GetSeckeyAsElement(sk)
	fmt.Println("recovery Ok? ", ckks.IsSame(ptxt1.El(), skEle))
}

func printDebugBrief(params *ckks.Parameters, valuesWant []complex128, valuesTest []complex128) {

	precStats := ckks.GetPrecisionStats(params, nil, nil, valuesWant, valuesTest)

	fmt.Println(precStats.String())
}

func printDebug(params *ckks.Parameters, ciphertext *ckks.Ciphertext, valuesWant []complex128, decryptor ckks.Decryptor, encoder ckks.Encoder) (valuesTest []complex128) {

	valuesTest = encoder.Decode(decryptor.DecryptNew(ciphertext), params.LogSlots())

	fmt.Println()
	fmt.Printf("Level: %d (logQ = %d)\n", ciphertext.Level(), params.LogQLvl(ciphertext.Level()))
	fmt.Printf("Scale: 2^%f\n", math.Log2(ciphertext.Scale()))
	fmt.Printf("ValuesTest: %6.10f %6.10f %6.10f %6.10f...\n", valuesTest[0], valuesTest[1], valuesTest[2], valuesTest[3])
	fmt.Printf("ValuesWant: %6.10f %6.10f %6.10f %6.10f...\n", valuesWant[0], valuesWant[1], valuesWant[2], valuesWant[3])
	fmt.Println()

	precStats := ckks.GetPrecisionStats(params, nil, nil, valuesWant, valuesTest)

	fmt.Println(precStats.String())

	return
}

func main() {
	example()
}
