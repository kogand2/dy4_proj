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
#include "mono_path.h"

// testing time complexity
#include <chrono>

void process_block_stream(int mode){
  auto start_time = std::chrono::high_resolution_clock::now();

  // Determine custom parameters depending on mode
  float rf_Fs;
  int rf_decim, audio_decim, audio_exp;

  switch(mode) {
    case 0:
      rf_Fs = 2400000.0;
      rf_decim = 10;
      audio_exp = 1;
      audio_decim = 5;
      break;
    case 1:
      rf_Fs = 1440000.0;
      rf_decim = 5;
      audio_exp = 1;
      audio_decim = 6;
      break;
    case 2:
      rf_Fs = 2400000.0;
      rf_decim = 10;
      audio_exp = 147;
      audio_decim = 800;
      break;
    case 3:
      rf_Fs = 2304000.0;
      rf_decim = 9;
      audio_exp = 441;
      audio_decim = 2560;
      break;
  }

  // RF front end variables
	float rf_Fc = 100000.0;
	int rf_taps = 151;

  // audio path variables
	float audio_Fc = 16000;
	int audio_taps = 151;

  // RF LPF filter coefficients
	std::vector<float> rf_coeff;
	low_pass_coeff(rf_Fs, rf_Fc, rf_taps, rf_coeff);

  // audio LPF filter coefficients
	std::vector<float> audio_coeff;
	low_pass_coeff((rf_Fs/rf_decim)*audio_exp, audio_Fc, audio_taps*audio_exp, audio_coeff);

	float block_size = 102400;

  // filtered data
  std::vector<float> i_filt;
	std::vector<float> q_filt;

  // state saving variable for audio data convolution
	std::vector<float> audio_state;
  audio_state.resize(audio_coeff.size() - 1, 0.0);

  // state saving variables for I and Q samples convolution
	std::vector<float> state_i_lpf_100k;
	std::vector<float> state_q_lpf_100k;
	state_i_lpf_100k.resize(rf_coeff.size() - 1, 0.0);
	state_q_lpf_100k.resize(rf_coeff.size() - 1, 0.0);

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

      // timing analysis
      auto stop_time = std::chrono::high_resolution_clock::now();
		  std::chrono::duration<double, std::milli> DFT_run_time = stop_time-start_time;
		  std::cerr << DFT_run_time.count() << " milliseconds" << "\n";

      std::cerr << "Use the following command to play the audio: \n";
      if (mode == 2 || mode == 3)
        std::cerr << "cat test.bin | aplay -c 1 -f S16_LE -r 44100 \n";
      else
        std::cerr << "cat test.bin | aplay -c 1 -f S16_LE -r 48000 \n";

      exit(1);
    }

    //std::cerr << "Processing block " << block_id << std::endl;

		// STEP 1: IQ samples demodulation
		std::vector<float> i_data(block_size / 2);
    std::vector<float> q_data(block_size / 2);

		for (int k = 0; k < block_size / 2; k++) {
      i_data[k] = block_data[2*k];
      q_data[k] = block_data[2*k + 1];
    }

    // filter out IQ data with convolution
		ds_block_conv(i_filt,i_data,rf_coeff,state_i_lpf_100k,rf_decim);
		ds_block_conv(q_filt,q_data,rf_coeff,state_q_lpf_100k,rf_decim);

		// take downsampled filtered IQ data
		std::vector<float>i_ds;
    i_ds.resize(i_filt.size()/rf_decim);
		std::vector<float>q_ds;
    q_ds.resize(q_filt.size()/rf_decim);

		downsample(rf_decim,i_filt,i_ds);
		downsample(rf_decim,q_filt,q_ds);

		// perform demodulation on IQ data
		IQ_demod = fmDemod(i_ds, q_ds, demod_state);
    //std::cerr << "test1\n";
		// STEP 2: Mono path
    std::vector<float> audio_block = mono_path(mode, IQ_demod, audio_coeff, audio_state, audio_decim, audio_exp);

    // STEP 3: prepare audio data for output
    std::vector<short int> audio_data;

    for (unsigned int k = 0; k < audio_block.size(); k++) {
      if(std::isnan(audio_block[k])){
        audio_data.push_back(0);
      }
      else{
        audio_data.push_back(static_cast<short int>(audio_block[k] * 16384));
      }
    }
    fwrite(&audio_data[0], sizeof(short int), audio_data.size(), stdout);
	}
}

int main(int argc, char* argv[])
{
  int mode = 0;

  // default operation
  if (argc < 2){
    std::cerr << "Operating in default mode 0" << std::endl;
  }
  // mode operation
  else if (argc == 2){
    mode = atoi(argv[1]);
    if (mode > 3){
      std::cerr << "Not a valid mode: " << mode << std::endl;
      exit(1);
    }
  }
  else{
    std::cerr << "Usage: " << argv[0] << std::endl;
    std::cerr << "or " << std::endl;
    std::cerr << "Usage: " << argv[0] << " <mode>" << std::endl;
    std::cerr << "\t\t <mode> is a value from 0 to 3" << std::endl;
    exit(1);
  }

  std::cerr << "Operating in mode " << mode << std::endl;
  process_block_stream(mode);

	return 0;
}
