# HEAAN
HEAAN defense is provided by Scheme::decryptForShare
which is not contained in any tagged version
the function that adds a gaussian noise to a polynomial is implemented in Ring::addGaussAndEqual
the function takes a sigma as parameter, and samples from a 2-D gaussian dist composed of two independent
gaussian distribution with stddev of sigma
the defense function(Scheme.cpp line 333 Scheme::decryptForShare) is described as below:
1. decrypt message
2. sigma1 = 3.2*sqrt(2)
3. if logErrBound == -1:
	sigma2 = $8*sigma1 / \sqrt{2pi*sigma1^2-64}\approx 4.5$
	else:
	sigma2 = $2^{logErrBound} / \sqrt{2pi}$
	(note that logErrBound defaults to 0)
4. add gaussian noise with stddev of sigma2 to each coefficient of message

# lattigo
the defense is implemented in branch dev_indCPA+\_migration, which is not merged into the main branch yet. so we are not going to discuss it

# HElib
1. an upper bound on ciphertext error is recorded in the library, let it be $\varepsilon$
2. the precision set by user is denoted by $\varepsilon'=2^{-k}$
3. if $\varepsilon'<\varepsilon$, a warning is raised, but anyway, the user-defined precision is used
4. $\sigma_{min}$ = 3.2\*2
5. $B=\sqrt{n\cdot lnn}$
6. we want $B\cdot \sigma_t / f = \varepsilon$, so the target stddev is $\varepsilon\cdot f / B$
7. if $\sigma_t < \sigma_{min}$, a warning is raised, the greater one among the two stddevs is used
8. seed(sk, c)
9. $e\leftarrow N(0,\sigma^2)$, then perform canonical embedding. rejection sampling for 1000 times
# PALISADE
Decoding:
defense is already present
the stddev in Gaussian Distribution is sqrt(2) * (the stddev of the image part in message vector)
the image part is cleared at the end of decoding for defense

NOTE: defense noise
defense noise is added in the following manner
suppose decoded vector z corresponds to the ring element msg = (m_0, m_1, ..., m_(n-1))
then conj(z) corresponds to msg' = (m_0, -m_(n-1), ..., -m_1)
2 * Image(z) = z - conj(z) corresponds to msg - msg' = (0, m_1 + m_(n-1), m_2 + m_(n-2), ..., m_1 + m_(n-1))
since only REAL vectors are accepted as inputs, the image part of z can somehow represent the error
let the stddev be the standard deviation of the n-1 non-zero coefficients in (msg - msg')/2
then the standard deviation for gaussian sampling is 2 * stddev
besides, 2 * Real(z) = z + conj(z) corresponds to msg + msg', which is used as the actual input for decoding
first (msg + msg')/2 is computed, then gaussian noise is added to each coefficient independently
Full description of the defense
(the polynomial coefficients are already in multi-precision representation)
1. the coefficients are rescaled s.t. the rescaling factor is exactly 2^p (where p is given by the user)
2. compute the stddev of the n-1 non-zero coefficients of (msg - msg')/2
3. if stddev < sqrt(n)/8, stddev := sqrt(n)/8 FIXME: why?
4. if log2(stddev) > p - 10, throw an exception that precision is less than 10 bits FIXME: why?
5. stddev := sqrt(M + 1) * stddev, where M defaults to 1
6. compute (msg + msg')/2
7. add gaussian noise to each coefficient in the polynomial above
8. perform decoding fft (canonical embedding)
9. clear the image part of the decoded vector
10. estimate the approximate error as round(log2(stddev * sqrt(n)))
