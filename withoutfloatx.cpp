#include "tfhe.h"
#include <iostream>

void Realfft(CArray& x, CArray& y)
{
	const size_t N = x.size();
	if (N <= 1) return;

	// divide
	CArray even = x[std::slice(0, N / 2, 2)];
	CArray  odd = x[std::slice(1, N / 2, 2)];

	// conquer
	Realfft(even, even);
	Realfft(odd, odd);

	// combine
	for (size_t k = 0; k < N / 2; ++k)
	{
		Complex a = Complex(std::polar(1.0, -2 * PI * k / N)) * odd[k];
		y[k] = even[k] + a;
		y[k + N / 2] = even[k] - a;
	}
}

void fft(TorusPolynomial* A, const LagrangeHalfCPolynomial* B) {
	Complex* test = new Complex[A->N];
	for (int i = 0; i < A->N; i++)
		test[i] = *((Complex*)B->data + i);
	CArray data(test, A->N);
	CArray result;
	result.resize(A->N);
	Realfft(data, result);
	Complex k(0.5, 0);
	for (int i = 0; i < A->N; i++) {
		if (real(result[i]) < 0)
			* (A->coefsT + i) = int32_t(real(result[i] - k));
		else
			*(A->coefsT + i) = int32_t(real(result[i] + k));
	}
}

// inverse fft (in-place)
void Realifft(int n, CArray& x, CArray& y) {
	{
		std::cout << n << std::endl;
		// n=1 return
		if (n == 1) {
			y[0] = x[0];
			return;
		}
		CArray A, B, C, D;
		A.resize(n / 2);
		B.resize(n / 2);
		C.resize(n / 2);
		D.resize(n / 2);

		//divide
		for (int k = 0; k < n / 2; ++k) {
			Complex a(0.5, 0);
			Complex b = Complex(std::polar(0.5, 2 * PI * k / n));
			A[k] = x[k] * a + x[k + n / 2] * a;
			B[k] = (b * x[k] - x[k + n / 2] * b);
		}
		//couquer
		Realifft(n / 2, A, C);
		Realifft(n / 2, B, D);

		//combine
		for (int k = 0; k < n / 2; ++k) {
			y[2 * k] = C[k];
			y[2 * k + 1] = D[k];
		}
	}
}
void ifft(int n, LagrangeHalfCPolynomial* a, IntPolynomial* b) {
	Complex* test = new Complex[n];
	for (int i = 0; i < n; i++)
		test[i] = Complex(*(b->coefs + i));
	CArray data(test, n);
	CArray result;
	result.resize(n);
	Realifft(n, data, result);
	for (int i = 0; i < n; i++)
		* ((Complex*)a->data + i) = result[i];
}

void tLweToFFTConvert(TLweSampleFFT* result, const TLweSample* source, const TLweParams* params) {
	const int32_t k = params->k;

	for (int32_t i = 0; i <= k; ++i)
		fft(result->a + i, source->a + i);
	result->current_variance = source->current_variance;
}

void tGswFFTExternMulToTLwe(TLweSample* accum, const TGswSampleFFT* gsw, const TGswParams* params) {
	const TLweParams* tlwe_params = params->tlwe_params;
	const int32_t k = tlwe_params->k;
	const int32_t l = params->l;
	const int32_t kpl = params->kpl;
	const int32_t N = tlwe_params->N;
	//TODO attention, improve these new/delete...
	IntPolynomial* deca = new_IntPolynomial_array(kpl, N); //decomposed accumulator
	LagrangeHalfCPolynomial* decaFFT = new_LagrangeHalfCPolynomial_array(kpl, N); //fft version
	TLweSampleFFT* tmpa = new_TLweSampleFFT(tlwe_params);

	for (int32_t i = 0; i <= k; i++)
		tGswTorus32PolynomialDecompH(deca + i * l, accum->a + i, params);
	for (int32_t p = 0; p < kpl; p++)
		ifft(N, decaFFT + p, deca + p);

	tLweFFTClear(tmpa, tlwe_params);
	for (int32_t p = 0; p < kpl; p++) {
		tLweFFTAddMulRTo(tmpa, decaFFT + p, gsw->all_samples + p, tlwe_params);
	}
	tLweFromFFTConvert(accum, tmpa, tlwe_params);

	delete_TLweSampleFFT(tmpa);
	delete_LagrangeHalfCPolynomial_array(kpl, decaFFT);
	delete_IntPolynomial_array(kpl, deca);
}



int main() {
	const int32_t nb_samples = 50;
	const Torus32 mu_boot = modSwitchToTorus32(1, 8);

	// generate params 
	int32_t minimum_lambda = 100;
	TFheGateBootstrappingParameterSet* params = new_default_gate_bootstrapping_parameters(minimum_lambda);
	const LweParams* in_out_params = params->in_out_params;
	// generate the secret keyset
	TFheGateBootstrappingSecretKeySet* keyset = new_random_gate_bootstrapping_secret_keyset(params);

	// generate input samples
	LweSample* test_in = new_LweSample_array(nb_samples, in_out_params);
	for (int32_t i = 0; i < nb_samples; ++i) {
		lweSymEncrypt(test_in + i, modSwitchToTorus32(i, nb_samples), 0.01, keyset->lwe_key);
	}
	// output samples
	LweSample* test_out = new_LweSample_array(nb_samples, in_out_params);


	// bootstrap input samples
	for (int32_t i = 0; i < nb_samples; ++i) {
		tfhe_bootstrap_FFT(test_out + i, keyset->cloud.bkFFT, mu_boot, test_in + i);
	}

	return 0;
}