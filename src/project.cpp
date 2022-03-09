/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include "block_conv_fn.h"


void downsample(int D,std::vector<float> input, std::vector<float> &down)
{
    for(unsigned int i = 0; i < input.size(); i = i + D)
    {
        down.push_back(input[i]);
    }
}

int main()
{
  std::cout << "test\n";
	// binary files can be generated through the
	// Python models from the "../model/" sub-folder
	const std::string in_fname = "../data/samples0.raw";
	std::vector<float> bin_data;
	readBinData(in_fname, bin_data);

	//change these with if statements based on modes
	//currently mode 0
	float rf_Fs = 2400000.0, rf_Fc = 100000.0;
	int rf_taps = 151, rf_decim = 10, rf_up = 0;

	float audio_Fs = 48000, audio_Fc = 16000;
	int audio_taps = 151, audio_decim = 5, audio_up = 0;

	//IQ data from -1 to +1 in 32-bit format
	std::vector<float> iq_data;
	iq_data.resize(bin_data.size(), 0.0);
	for(unsigned int i = 0; i < bin_data.size(); i++)
	{
			iq_data[i] = (bin_data[i] - 128.0)/128.0;
	}
  printRealVector(iq_data);
	std::vector<float> rf_coeff;
	low_pass_coeff(rf_Fs, rf_Fc, rf_taps, rf_coeff);

	std::vector<float> audio_coeff;
	low_pass_coeff(rf_Fs/rf_decim, audio_Fc, audio_taps, audio_coeff);

	float block_size = 1024*rf_decim*audio_decim*2;
	int block_count = 0;


  std::vector<float> i_filt;
	std::vector<float> q_filt;
  std::vector<float> audio_filt;
	std::vector<float> state_i_lpf_100k;
	std::vector<float> state_q_lpf_100k;
	state_i_lpf_100k.resize(rf_taps - 1, 0.0);
	state_q_lpf_100k.resize(rf_taps - 1, 0.0);

	std::vector<float> dummy_state;
	std::vector<float> dummy_fm;
	std::vector<float> filt_block;
	std::vector<float> audio_data;

	filt_block.resize(audio_coeff.size(), 0.0);
	dummy_state.resize(2, 0.0);


	while ((block_count+1)*block_size < iq_data.size())
	{
		// printf('Processing block %d\n',block_count);
		// i and q block convolution
		std::vector<float>block_split_i;
    std::vector<float>block_split_q;
		std::vector<float>block_split_down_i;
    std::vector<float>block_split_down_q;


		block_split_i = std::vector<float>(iq_data.begin() + block_count*block_size, iq_data.begin() + (block_count + 1)*block_size);
    block_split_q = std::vector<float>(iq_data.begin() + block_count*block_size + 1, iq_data.begin() + (block_count + 1)*block_size);

		downsample(2, block_split_i, block_split_down_i);
    downsample(2, block_split_q, block_split_down_q);

		state_block_conv(i_filt,block_split_down_i,rf_coeff,state_i_lpf_100k);

		state_block_conv(q_filt,block_split_down_q,rf_coeff,state_q_lpf_100k);

		//decimation
		std::vector<float>i_ds;
		std::vector<float>q_ds;

		downsample(rf_decim,i_filt,i_ds);
		downsample(rf_decim,q_filt,q_ds);

		//demodulation
		dummy_fm = fmDemod(i_ds, q_ds, dummy_state);

		//block convolution
		state_block_conv(audio_filt, dummy_fm, audio_coeff, filt_block);

		std::vector<float>audio_block;
		downsample(audio_decim,audio_filt,audio_block);
    std::cout << audio_filt.size() << "\n";
    std::cout << audio_block.size() << "\n";

		audio_data.insert(audio_data.end(), audio_block.begin(), audio_block.end());
		block_count += 1;
	}
  //printRealVector(audio_data);
	const std::string out_fname = "../data/float32filtered.bin";
	writeBinData(out_fname, audio_data);

	// generate an index vector to be used by logVector on the X axis
	std::vector<float> vector_index;
	genIndexVector(vector_index, bin_data.size());
	// log time data in the "../data/" subfolder in a file with the following name
	// note: .dat suffix will be added to the log file in the logVector function
	logVector("demod_time", vector_index, bin_data);

	// take a slice of data with a limited number of samples for the Fourier transform
	// note: NFFT constant is actually just the number of points for the
	// Fourier transform - there is no FFT implementation ... yet
	// unless you wish to wait for a very long time, keep NFFT at 1024 or below
	std::vector<float> slice_data = \
		std::vector<float>(bin_data.begin(), bin_data.begin() + NFFT);
	// note: make sure that binary data vector is big enough to take the slice

	// declare a vector of complex values for DFT
	std::vector<std::complex<float>> Xf;
	// ... in-lab ...
	// compute the Fourier transform
	// the function is already provided in fourier.cpp

	// compute the magnitude of each frequency bin
	// note: we are concerned only with the magnitude of the frequency bin
	// (there is NO logging of the phase response)
	std::vector<float> Xmag;
	// ... in-lab ...
	// compute the magnitude of each frequency bin
	// the function is already provided in fourier.cpp

	// log the frequency magnitude vector
	vector_index.clear();
	genIndexVector(vector_index, Xmag.size());
	logVector("demod_freq", vector_index, Xmag); // log only positive freq

	// for your take-home exercise - repeat the above after implementing
	// your OWN function for PSD based on the Python code that has been provided
	// note the estimate PSD function should use the entire block of "bin_data"
	//
	// ... complete as part of the take-home ...
	//

	// if you wish to write some binary files, see below example
	//
	// const std::string out_fname = "../data/outdata.bin";
	// writeBinData(out_fname, bin_data);
	//
	// output files can be imported, for example, in Python
	// for additional analysis or alternative forms of visualization

	// naturally, you can comment the line below once you are comfortable to run GNU plot
	std::cout << "Run: gnuplot -e 'set terminal png size 1024,768' ../data/example.gnuplot > ../data/example.png\n";

	return 0;
}
