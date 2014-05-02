#include <stdio.h>
#include <unistd.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <complex>
#include <valarray>
#include <plplot/plstream.h>

#define L 10
#define dx 0.1
#define dt 0.000001

#define h 1
#define m 1

const int n=L/dx;
typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;
plstream *pls;

CArray psi(n);
double x[n]={};
double V[n] = {};

const double pi = 3.141592653589793238460;

void fft(CArray& x);
void ifft(CArray& x);
void initpsi();
void plot();
void setpotential();
void normrescale();
void applypotential();
void applymomentum();

using namespace std;

int main(int argc, char * argv[])
{
	cout << endl << "Length: " << L << endl;
	cout << "dx: " << dx << endl;
	cout << "N: " << n << endl << endl;
	
	for (int i=0;i<n;i++)
		x[i] = dx * i;
	initpsi();
	
	pls = new plstream();
	pls -> sdev("xcairo");
	pls -> init();
	pls -> env(0,L,0,3,0,1);

	//Setting norm to unity
	normrescale();

	//Creating potential
	setpotential();
	
	//Iterative loop
	double normsum;
	for (int i=0; i<10000; i++)
	{
	plot();
	//applypotential();
	fft(psi);
	normsum = 0;
	for (int i = 0;i<n;i++)
		normsum += norm(psi[i]) * dx ;
	cout << "norm of wavefunction=" << normsum << endl;
	sleep(1);

	//applymomentum();
	ifft(psi);
	normsum = 0;
	for (int i = 0;i<n;i++)
		normsum += norm(psi[i]) * dx ;
	cout << "norm of wavefunction=" << normsum << endl;
	sleep(1);
	}
	//plot();
	return 0;
}

void applymomentum()
{
	for (int i=0;i<n;i++)
		psi[i] *= exp(Complex(0,-dt*i*i/(2*m*h)));
}

void applypotential()
{
	for (int i=0;i<n;i++)
	{
		psi[i] *= exp(Complex(0,-dt*V[i]/h));
	}
}

void setpotential()
{


}

void normrescale()
{
	double normsum = 0;
	for (int i = 0;i<n;i++)
		normsum += norm(psi[i]) * dx ;
	cout << "norm of wavefunction=" << normsum << "\tNow correcting" << endl;
	for (int i = 0;i<n;i++)
		psi[i]/=sqrt(normsum);
}

void plot()
{
	double psinorm[n];
	for (int i=0;i<n;i++)
		psinorm[i] = abs(psi[i]);
	pls -> clear();
	pls -> line(n,x,psinorm);
	pls ->flush();
}	

void initpsi()
{
	/*double b=sqrt(3)/4.;
	  for (int i=20;i<40;i++)
	  {
	  psi[i] = (i-20)* b * dx;
	  }
	  for (int i=40;i<60;i++)
	  {
	  psi[i] = (60-i)* b  * dx;
	  }*/
	for (int i=0; i<n;i++)
	{
		psi[i] = exp(-pow((dx*i-L/4),2)/2);
	}

}


// Cooley–Tukey FFT (in-place)
void fft(CArray& x)
{
	const size_t N = x.size();
	if (N <= 1) return;

	// divide
	CArray even = x[std::slice(0, N/2, 2)];
	CArray  odd = x[std::slice(1, N/2, 2)];

	// conquer
	fft(even);
	fft(odd);

	// combine
	for (size_t k = 0; k < N/2; ++k)
	{
		Complex t = std::polar(1.0, -2 * pi * k / N) * odd[k];
		x[k    ] = even[k] + t;
		x[k+N/2] = even[k] - t;
	}
}

// inverse fft (in-place)
void ifft(CArray& x)
{
	// conjugate the complex numbers
	x = x.apply(std::conj);

	// forward fft
	fft( x );

	// conjugate the complex numbers again
	x = x.apply(std::conj);

	// scale the numbers
	x /= x.size();
}
