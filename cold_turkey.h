
#ifndef COLD_TURKEY_INCLUDED
#define COLD_TURKEY_INCLUDED


#include <stdint.h>



/*
calculates DFT of signal and stores it in freqs
may be faster than fft for short signals

signal_real: real component of input signal
signal_imaginary: imaginary component of input signal (usually all 0s)
freqs_real: real component of output frequencies
freqs_imaginary: imaginary component of output frequencies
length: the length of all the arrays (they must all be the same length)
*/
void dft(double* signal_real, double* signal_imaginary, double* freqs_real, double* freqs_imaginary, uint32_t length);

/*
calculates inverse DFT of freqs and stores it in signal
may be faster than ifft for short inputs

freqs_real: real component of input frequencies
freqs_imaginary: imaginary component of input frequencies
signal_real: real component of output signal
signal_imaginary: imaginary component of output signal
length: the length of all the arrays (they must all be the same length)
*/
void idft(double* freqs_real, double* freqs_imaginary, double* signal_real, double* signal_imaginary, uint32_t length);


/*
an in-place FFT
length must be a power of 2 and less than 2^32
returns -1 if unsuccessful

signal_real: real component of input signal
signal_imaginary: imaginary component of input signal (usually all 0s)
length: the length of the signal
*/
int fft(double* signal_real, double* signal_imaginary, uint32_t length);

/*
an in-place inverse FFT
length must be a power of 2 and less than 2^32
returns -1 if unsuccessful

freqs_real: real component of input frequencies
freqs_imaginary: imaginary component of input frequencies
length: the length of the frequency arrays
*/
int ifft(double* freqs_real, double* freqs_imaginary, uint32_t length);


/*
calculates FFT of signal and stores it in freqs
length must be a power of 2 and less than 2^32
returns -1 if unsuccessful

signal_real: real component of input signal
signal_imaginary: imaginary component of input signal (usually all 0s)
freqs_real: real component of output frequencies
freqs_imaginary: imaginary component of output frequencies
length: the length of all the arrays (they must all be the same length)
*/
int fft_out_of_place(double* signal_real, double* signal_imaginary, double* freqs_real, double* freqs_imaginary, uint32_t length);

/*
calculates inverse FFT of freqs and stores it in signal
length must be a power of 2 and less than 2^32
returns -1 if unsuccessful

freqs_real: real component of input frequencies
freqs_imaginary: imaginary component of input frequencies
signal_real: real component of output signal
signal_imaginary: imaginary component of output signal
length: the length of all the arrays (they must all be the same length)
*/
int ifft_out_of_place(double* freqs_real, double* freqs_imaginary, double* signal_real, double* signal_imaginary, uint32_t length);




/*
------------------------- start of float function definitions --------------------------
*/


/*
calculates DFT of signal and stores it in freqs
may be faster than fft for short signals

signal_real: real component of input signal
signal_imaginary: imaginary component of input signal (usually all 0s)
freqs_real: real component of output frequencies
freqs_imaginary: imaginary component of output frequencies
length: the length of all the arrays (they must all be the same length)
*/
void dft_float(float* signal_real, float* signal_imaginary, float* freqs_real, float* freqs_imaginary, uint32_t length);

/*
calculates inverse DFT of freqs and stores it in signal
may be faster than ifft for short inputs

freqs_real: real component of input frequencies
freqs_imaginary: imaginary component of input frequencies
signal_real: real component of output signal
signal_imaginary: imaginary component of output signal
length: the length of all the arrays (they must all be the same length)
*/
void idft_float(float* freqs_real, float* freqs_imaginary, float* signal_real, float* signal_imaginary, uint32_t length);


/*
an in-place FFT
length must be a power of 2 and less than 2^32
returns -1 if unsuccessful

signal_real: real component of input signal
signal_imaginary: imaginary component of input signal (usually all 0s)
length: the length of the signal
*/
int fft_float(float* signal_real, float* signal_imaginary, uint32_t length);

/*
an in-place inverse FFT
length must be a power of 2 and less than 2^32
returns -1 if unsuccessful

freqs_real: real component of input frequencies
freqs_imaginary: imaginary component of input frequencies
length: the length of the frequency arrays
*/
int ifft_float(float* freqs_real, float* freqs_imaginary, uint32_t length);


/*
calculates FFT of signal and stores it in freqs
length must be a power of 2 and less than 2^32
returns -1 if unsuccessful

signal_real: real component of input signal
signal_imaginary: imaginary component of input signal (usually all 0s)
freqs_real: real component of output frequencies
freqs_imaginary: imaginary component of output frequencies
length: the length of all the arrays (they must all be the same length)
*/
int fft_out_of_place_float(float* signal_real, float* signal_imaginary, float* freqs_real, float* freqs_imaginary, uint32_t length);

/*
calculates inverse FFT of freqs and stores it in signal
length must be a power of 2 and less than 2^32
returns -1 if unsuccessful

freqs_real: real component of input frequencies
freqs_imaginary: imaginary component of input frequencies
signal_real: real component of output signal
signal_imaginary: imaginary component of output signal
length: the length of all the arrays (they must all be the same length)
*/
int ifft_out_of_place_float(float* freqs_real, float* freqs_imaginary, float* signal_real, float* signal_imaginary, uint32_t length);



#endif