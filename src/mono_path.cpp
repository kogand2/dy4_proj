#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include "block_conv_fn.h"
#include "mono_path.h"

// read in raw block data
void readStdinBlockData(unsigned int num_samples, std::vector<float> &block_data){
  std::vector<char> raw_data(num_samples);
  std::cin.read(reinterpret_cast<char*>(&raw_data[0]), num_samples*sizeof(char));
  for(int k = 0; k < (int)num_samples; k++) {
    block_data[k] = float(((unsigned char)raw_data[k] - 128)/ 128.0);
  }
}

// dowsampling function
void downsample(int D,std::vector<float> input, std::vector<float> &down)
{
    for(unsigned int i = 0; i < input.size(); i = i + D)
    {
        down.push_back(input[i]);
    }
}

// mono path
void mono_path(int mode)
{
  // if mode == ...
  // RF front end variables
	float rf_Fs = 2400000.0, rf_Fc = 100000.0;
	int rf_taps = 151, rf_decim = 10, rf_up = 0;

  // audio path variables
	float audio_Fs = 48000, audio_Fc = 16000;
	int audio_taps = 151, audio_decim = 5, audio_up = 0;

  // RF LPF filter coefficients
	std::vector<float> rf_coeff;
	low_pass_coeff(rf_Fs, rf_Fc, rf_taps, rf_coeff);

  // audio path LPF coefficients
	std::vector<float> audio_coeff;
	low_pass_coeff(rf_Fs/rf_decim, audio_Fc, audio_taps, audio_coeff);

	float block_size = 1024*rf_decim*audio_decim*2;
	int block_count = 0;

  // filtered data
  std::vector<float> i_filt;
	std::vector<float> q_filt;
  std::vector<float> audio_filt;

  // state saving variable for audio data convolution
	std::vector<float> audio_state;
  audio_state.resize(audio_coeff.size(), 0.0);

  // state saving variables for I and Q samples convolution
	std::vector<float> state_i_lpf_100k;
	std::vector<float> state_q_lpf_100k;
	state_i_lpf_100k.resize(rf_taps - 1, 0.0);
	state_q_lpf_100k.resize(rf_taps - 1, 0.0);

  // demodulation variables
	std::vector<float> demod_state;
	std::vector<float> IQ_demod;
	demod_state.resize(2, 0.0);

  // decipher each block
	for(unsigned int block_id = 0; ; block_id++) {
		std::vector<float> block_data(block_size);
    readStdinBlockData(block_size, block_data);
    if ((std::cin.rdstate()) != 0) {
      std::cerr << "End of input stream reached" << std::endl;
      exit(1);
    }
    std::cerr << "Read block " << block_id << std::endl;

		// STEP 1: IQ samples demodulation
		std::vector<float> i_data(block_size / 2);
    std::vector<float> q_data(block_size / 2);

		for (int k = 0; k < block_size / 2; k++) {
      i_data[k] = block_data[2*k];
      q_data[k] = block_data[2*k + 1];
    }

    // filter out IQ data with convolution
		state_block_conv(i_filt,i_data,rf_coeff,state_i_lpf_100k);
		state_block_conv(q_filt,q_data,rf_coeff,state_q_lpf_100k);

		// downsample filtered IQ data
		std::vector<float>i_ds;
		std::vector<float>q_ds;

		downsample(rf_decim,i_filt,i_ds);
		downsample(rf_decim,q_filt,q_ds);

		// perform demodulation on IQ data
		IQ_demod = fmDemod(i_ds, q_ds, demod_state);

		// STEP 2: Mono path

    // filter out audio data with convolution
		state_block_conv(audio_filt, IQ_demod, audio_coeff, audio_state);

    // downsample filtered audio data
		std::vector<float>audio_block;
		downsample(audio_decim,audio_filt,audio_block);

    // STEP 3: prepare audio data for output
    std::vector<short int> audio_data;

    for (unsigned int k = 0; k < audio_block.size(); k++) {
      if(std::isnan(audio_block[k])){
        audio_data.push_back(0);
      }else{
        audio_data.push_back(static_cast<short int>(audio_block[k] * 16384));
      }
    }
    fwrite(&audio_data[0], sizeof(short int), audio_data.size(), stdout);
	}
}
