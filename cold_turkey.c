#include <stdint.h>
#include <math.h>

#define PI 3.14159265358979323846



// multiplies a by b and returns the real component
#define complex_multiply_real(a_real, a_imaginary, b_real, b_imaginary) \
    ((a_real) * (b_real) - (a_imaginary) * (b_imaginary))

// multiplies a by b and returns the imaginary component
#define complex_multiply_imaginary(a_real, a_imaginary, b_real, b_imaginary) \
    ((a_real) * (b_imaginary) + (a_imaginary) * (b_real))



void dft_loops(double* signal_real, double* signal_imaginary, double* freqs_real, double* freqs_imaginary, 
    uint32_t length, int direction, uint32_t normalization_factor) 
{
    for (uint32_t f = 0; f < length; f++) {
        double real_sum = 0;
        double imaginary_sum = 0;
        for (uint32_t t = 0; t < length; t++) {
            
            // evaluate real and imaginary components of e^(-2*i*pi*f*t/length) individually using euler's formula
            double exponent = direction * 2 * PI * f * t / length;
            double exponent_real = cos(exponent);
            double exponent_imaginary = sin(exponent);

            // binomial expansion for multiplying two complex numbers
            real_sum += complex_multiply_real(signal_real[t], signal_imaginary[t], exponent_real, exponent_imaginary);
            imaginary_sum += complex_multiply_imaginary(signal_real[t], signal_imaginary[t], exponent_real, exponent_imaginary);
        }
        freqs_real[f] = real_sum / normalization_factor;
        freqs_imaginary[f] = imaginary_sum / normalization_factor;
    }
}


// edited version of bit reverse function from the algorithm-archive: https://github.com/algorithm-archivists/algorithm-archive/blob/master/contents/cooley_tukey/code/c%2B%2B/fft.cpp
void sort_in_bit_reversed_order(double* array, uint32_t length) {
    uint32_t num_bits = (uint32_t)(log2(length)); // number of bits needed to represent any location in array
    uint32_t num_swaps = (uint32_t)(ceil(log2(num_bits))); // this algorithm works by repeatedly shifting around smaller and smaller groups of bits, we COULD do all the bit group swaps for all input sizes but it's a waste of time if we only need to reverse like 5 bits
    uint32_t num_trailing_zeros = (uint32_t)(pow(2, num_swaps) - num_bits); // number of trailing 0's that we need to get rid of when the bit string is reversed

    for (int index = 0; index < length; index++) {

        uint32_t rev_index = index;                                                                 // 00000012 3456789A (example 10 bit index)
        switch (num_swaps) { // switch is just for efficiency: do as few bit swaps as possible
            case 5: // start here if num_bits > 16
                rev_index = ((rev_index >> 16) | (rev_index << 16));                                // swap neighboring groups of 16 bits (for above 16 bit input only so skipped in 10 bit example)
            case 4: // start here if num_bits > 8
                rev_index = (((rev_index & 0xff00ff00) >> 8) | ((rev_index & 0x00ff00ff) << 8));    // 3456789A 00000012 (swap neighboring groups of 8 bits)
            case 3: // start here if num_bits > 4
                rev_index = (((rev_index & 0xf0f0f0f0) >> 4) | ((rev_index & 0x0f0f0f0f) << 4));    // 789A3456 00120000 (swap neighboring groups of 4 bits)
            case 2: // start here if num_bits > 2
                rev_index = (((rev_index & 0xcccccccc) >> 2) | ((rev_index & 0x33333333) << 2));    // 9A785634 12000000 (swap neighboring groups of 2 bits)
            case 1: // start here if num_bits > 1
                rev_index = (((rev_index & 0xaaaaaaaa) >> 1) | ((rev_index & 0x55555555) << 1));    // A9876543 21000000 (swap neighboring bits)
        }
        rev_index = rev_index >> num_trailing_zeros;                                                // 000000A9 87654321 (shift to remove trailing 0s)
        
        if (rev_index > index) { // make sure we don't swap the same 2 elements twice
            double value_at_index = array[index];
            array[index] = array[rev_index];
            array[rev_index] = value_at_index;
        }
    }
}


// an iterative, decimation in time, radix 2 cooley tukey
int radix_2_cooley_tukey(double* freqs_real, double* freqs_imaginary, uint32_t total_length, int direction) {

    // throw error if direction or total_length is wrong
    if (((total_length & (total_length - 1)) != 0) || // check if total_length is a power of 2
        (direction != 1 && direction != -1)) // check if direction is clockwise/counterclockwise
    {
        return -1;
    }

    sort_in_bit_reversed_order(freqs_real, total_length);
    sort_in_bit_reversed_order(freqs_imaginary, total_length);

    // do length 2 DFTs first (seperate loop because skipping twiddle factors)
    for (uint32_t even_index = 0; even_index < total_length; even_index += 2) {
        uint32_t odd_index = even_index + 1;
        double even_value_real = freqs_real[even_index];
        double even_value_imaginary = freqs_imaginary[even_index];

        freqs_real[even_index] = even_value_real + freqs_real[odd_index];
        freqs_real[odd_index] = even_value_real - freqs_real[odd_index];
        freqs_imaginary[even_index] = even_value_imaginary + freqs_imaginary[odd_index];
        freqs_imaginary[odd_index] = even_value_imaginary - freqs_imaginary[odd_index];
    }
    
    // do all larger DFTs 
    for (uint32_t dft_length = 4; dft_length <= total_length; dft_length *= 2) {
        uint32_t num_dfts = total_length / dft_length;
        for (uint32_t dft_num = 0; dft_num < num_dfts; dft_num++) {
            uint32_t dft_start = dft_num * dft_length;
            uint32_t half_dft_length = dft_length / 2;
            for (uint32_t frequency_num = 0; frequency_num < half_dft_length; frequency_num++) {

                uint32_t even_index = dft_start + frequency_num;
                uint32_t odd_index = even_index + half_dft_length;
                double even_value_real = freqs_real[even_index];
                double even_value_imaginary = freqs_imaginary[even_index];
                
                // evaluate real and imaginary components of twiddle factor individually using euler's formula
                double exponent = direction * 2 * PI * frequency_num / dft_length;
                double twiddle_real = cos(exponent);
                double twiddle_imaginary = sin(exponent);

                // multiply twiddle factor by odd value
                double odd_times_twiddle_real = complex_multiply_real(
                    freqs_real[odd_index], freqs_imaginary[odd_index], twiddle_real, twiddle_imaginary);
                double odd_times_twiddle_imaginary = complex_multiply_imaginary(
                    freqs_real[odd_index], freqs_imaginary[odd_index], twiddle_real, twiddle_imaginary);

                freqs_real[even_index] = even_value_real + odd_times_twiddle_real;
                freqs_real[odd_index] = even_value_real - odd_times_twiddle_real;
                freqs_imaginary[even_index] = even_value_imaginary + odd_times_twiddle_imaginary;
                freqs_imaginary[odd_index] = even_value_imaginary - odd_times_twiddle_imaginary;
            }
        }
    }
    
    return 0;
}



void dft(double* signal_real, double* signal_imaginary, double* freqs_real, double* freqs_imaginary, uint32_t length) {
    dft_loops(signal_real, signal_imaginary, freqs_real, freqs_imaginary, length, -1, 1);
}


void idft(double* freqs_real, double* freqs_imaginary, double* signal_real, double* signal_imaginary, uint32_t length) {
    dft_loops(freqs_real, freqs_imaginary, signal_real, signal_imaginary, length, 1, length);
}


int fft(double* signal_real, double* signal_imaginary, uint32_t length) {
    return radix_2_cooley_tukey(signal_real, signal_imaginary, length, -1); // -1 for counterclockwise
}


int ifft(double* freqs_real, double* freqs_imaginary, uint32_t length) {
    int success = radix_2_cooley_tukey(freqs_real, freqs_imaginary, length, 1); // 1 for clockwise
    if (success != -1) {
        for (uint32_t sample = 0; sample < length; sample++) {
            freqs_real[sample] = freqs_real[sample] / length;
            freqs_imaginary[sample] = freqs_imaginary[sample] / length;
        }
    }
    return success;
}


int fft_out_of_place(double* signal_real, double* signal_imaginary, double* freqs_real, double* freqs_imaginary, uint32_t length) {
    // copy signal into frequency array and then do in place fft
    for (uint32_t sample = 0; sample < length; sample++) {
        freqs_real[sample] = signal_real[sample];
        freqs_imaginary[sample] = signal_imaginary[sample];
    }
    return fft(freqs_real, freqs_imaginary, length);
}


int ifft_out_of_place(double* freqs_real, double* freqs_imaginary, double* signal_real, double* signal_imaginary, uint32_t length) {
    // copy frequencies into signal array and then do in place ifft
    for (uint32_t bin = 0; bin < length; bin++) {
        signal_real[bin] = freqs_real[bin];
        signal_imaginary[bin] = freqs_imaginary[bin];
    }
    return ifft(signal_real, signal_imaginary, length);
}