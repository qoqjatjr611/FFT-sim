#include <complex>
#include <fstream>
#include <iostream>
#include <valarray>
#include "floatx.hpp"
#include <math.h>
#include <vector>
#include <random>

const double PI = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821;
typedef std::complex<flx::floatx<2, 36>> Complex;
typedef std::valarray<Complex> CArray;

// Cooley–Tukey FFT (in-place, divide-and-conquer)
// Higher memory requirements and redundancy although more intuitive
void fft(CArray& x)
{
	const size_t N = x.size();
	if (N <= 1) return;

	// divide
	CArray even = x[std::slice(0, N / 2, 2)];
	CArray  odd = x[std::slice(1, N / 2, 2)];

	// conquer
	fft(even);
	fft(odd);

	// combine
	for (size_t k = 0; k < N / 2; ++k)
	{
		Complex a = Complex(std::polar(1.0, -2 * PI * k / N)) * odd[k];
		x[k] = even[k] + a;
		x[k + N / 2] = even[k] - a;
	}
}

// inverse fft (in-place)
void ifft(CArray& x)
{
	// conjugate the complex numbers
	x = x.apply(std::conj);

	// forward fft
	fft(x);

	// conjugate the complex numbers again
	x = x.apply(std::conj);

	// scale the numbers
	for (int i = 0; i < x.size(); i++)
		x[i] /= x.size();
}
void fff(int points) {
	
	int count = 0;
	/*std::string filename("input.txt");
	std::ifstream input_file(filename);
	std::vector<Complex> nums;
	int number;


	while (input_file >> number) {
		nums.push_back(Complex(number));
	}

	Complex* test = new Complex[points];
	copy(nums.begin(), nums.end(), test);*/

	Complex* test = new Complex[points];
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<int> dis(0, 10000000000);
	for (int i = 0; i < points; i++) {
		test[i] = dis(gen) - 5000000000;
	}
	CArray data(test, points);
	/*
	for (int i = 0; i < points; ++i)
	{
		std::cout << data[i] << std::endl;
	}
	*/
	// forward fft
	fft(data);
	/*
	std::cout << "fft" << std::endl;
	for (int i = 0; i < points; ++i)
	{
		std::cout << data[i] << std::endl;
	}
	*/
	// inverse fft
	ifft(data);

	for (int i = 0; i < points; ++i)
	{
		if (fabs(double(real(data[i])) - double(real(test[i]))) >= 0.5 || fabs(double(imag(data[i])) - double(imag(test[i]))) >= 0.5)
			count++;
	}
	if (count != 0) {
		std::cout << "FFT FAIL!" << std::endl;
		exit(0);
	}
}
int main()
{
	int c;
	int points = 1;
	do {
		std::cout << "HOW MANY POINTS? (1~20): ";
		std::cin >> c;
	} while (c < 1 || c > 20);

	for (int i = 0; i < c; i++)
		points *= 2;
	for (int i = 0; i < 300; i++)
		fff(points);
	return 0;
}