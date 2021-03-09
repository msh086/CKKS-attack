package main

import (
	"fmt"
	"github.com/ldsec/lattigo/v2/ckks"
	"math"
	"math/cmplx"
	"time"
)

func example() {

	var start time.Time
	var err error

	LogN := uint64(14)
	LogSlots := uint64(13)

	LogModuli := ckks.LogModuli{
		LogQi: []uint64{55, 40, 40, 40, 40, 40, 40, 40},
		LogPi: []uint64{45, 45},
	}

	Scale := float64(1 << 40)

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

	r := float64(16)

	pi := 3.141592653589793

	slots := params.Slots()

	values := make([]complex128, slots)
	for i := range values {
		values[i] = complex(2*pi, 0)
	}


	/*
	   a plaintext polynomial is in NTT form, it's represented in a (nbModuli x N) matrix of uint64
	*/
	plaintext := ckks.NewPlaintext(params, params.MaxLevel(), params.Scale()/r)
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
		the modulus are deliberately chosen(not checked, just a guess) s.t. q_i = 1 mod N
		i.e. NTT modulus is the same as non-NTT modulus
		which means sampling a random element in NTT form is essentially the same as sampling a random element
		in non-NTT form
	*/
	ciphertext := encryptor.EncryptNew(plaintext)

	fmt.Printf("Done in %s \n", time.Since(start))

	printDebug(params, ciphertext, values, decryptor, encoder)

	fmt.Println()
	fmt.Println("===============================================")
	fmt.Printf("        EVALUATION OF i*x on %d values\n", slots)
	fmt.Println("===============================================")
	fmt.Println()

	start = time.Now()

	evaluator.MultByi(ciphertext, ciphertext)

	fmt.Printf("Done in %s \n", time.Since(start))

	for i := range values {
		values[i] *= complex(0, 1)
	}

	printDebug(params, ciphertext, values, decryptor, encoder)

	fmt.Println()
	fmt.Println("===============================================")
	fmt.Printf("       EVALUATION of x/r on %d values\n", slots)
	fmt.Println("===============================================")
	fmt.Println()

	start = time.Now()

	ciphertext.MulScale(r)

	fmt.Printf("Done in %s \n", time.Since(start))

	for i := range values {
		values[i] /= complex(r, 0)
	}

	printDebug(params, ciphertext, values, decryptor, encoder)

	fmt.Println()
	fmt.Println("===============================================")
	fmt.Printf("       EVALUATION of e^x on %d values\n", slots)
	fmt.Println("===============================================")
	fmt.Println()

	start = time.Now()

	coeffs := []complex128{
		complex(1.0, 0),
		complex(1.0, 0),
		complex(1.0/2, 0),
		complex(1.0/6, 0),
		complex(1.0/24, 0),
		complex(1.0/120, 0),
		complex(1.0/720, 0),
		complex(1.0/5040, 0),
	}

	poly := ckks.NewPoly(coeffs)

	if ciphertext, err = evaluator.EvaluatePoly(ciphertext, poly, rlk); err != nil {
		panic(err)
	}

	fmt.Printf("Done in %s \n", time.Since(start))

	for i := range values {
		values[i] = cmplx.Exp(values[i])
	}

	printDebug(params, ciphertext, values, decryptor, encoder)

	fmt.Println()
	fmt.Println("===============================================")
	fmt.Printf("       EVALUATION of x^r on %d values\n", slots)
	fmt.Println("===============================================")
	fmt.Println()

	start = time.Now()

	evaluator.Power(ciphertext, uint64(r), rlk, ciphertext)

	fmt.Printf("Done in %s \n", time.Since(start))

	for i := range values {
		values[i] = cmplx.Pow(values[i], complex(r, 0))
	}

	printDebug(params, ciphertext, values, decryptor, encoder)

	fmt.Println()
	fmt.Println("=========================================")
	fmt.Println("         DECRYPTION & DECODING           ")
	fmt.Println("=========================================")
	fmt.Println()

	start = time.Now()

	outPlain := decryptor.DecryptNew(ciphertext)
	outvalues := encoder.Decode(outPlain, params.LogSlots())

	fmt.Printf("Done in %s \n", time.Since(start))

	printDebug(params, ciphertext, values, decryptor, encoder)

	// now attack
	fmt.Println("log2(scale) = ", math.Log2(ciphertext.Scale()))
	fmt.Println("log2(inf norm) = ", )
	ptxt1 := ckks.NewPlaintext(params, ciphertext.Level(), ciphertext.Scale())
	ptxt2 := ckks.NewPlaintext(params, ciphertext.Level(), ciphertext.Scale())

	encoder.EncodeNTT(ptxt1, outvalues, params.LogSlots())
	fmt.Println("re-encoding ok? ", ckks.IsSame(ptxt1.El(), outPlain.El()))
	// now ptxt1 = m+e
	exposed := ckks.CastExposed(evaluator)
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
