#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include "block_conv_fn.h"
#include "mono_path.h"

#include <stdio.h>
#include <math.h>

// testing time complexity
#include <chrono>


void band_pass_coeff(float Fb, float Fe, float Fs, int num_taps, std::vector<float> &h)
{
  float normCenter = ((Fe + Fb)/2)/(Fs/2);
	float normPass = (Fe - Fb)/(Fs/2);

  h.clear();
	h.resize(num_taps, 0.0);

  for (int i = 0; i < num_taps; i++){
    if (i == (num_taps-1)/2)
      h[i] = normPass;
    else
      h[i] = normPass * (std::sin(PI*(normPass/2)*(i-((num_taps-1)/2))))/(PI*(normPass/2)*(i-((num_taps-1)/2)));
    h[i] = h[i] * std::cos(i*PI*normCenter);
    h[i] = 2.5 * h[i] * pow(std::sin(i*PI/num_taps), 2);
  }
}

void fmPll(std::vector<float> &pllIn, std::vector<float> &ncoOut, std::vector<float> &state, float freq, float Fs, float ncoScale = 2.0, float phaseAdjust = 0.0, float normBandwidth = 0.01)
{
    float Cp = 2.666;
    float Ci = 3.555;

    float Kp = normBandwidth*Cp;
    float Ki = normBandwidth*normBandwidth*Ci;

    ncoOut.clear();
  	ncoOut.resize(pllIn.size()+1);

    float integrator = state[0];
  	float phaseEst = state[1];
  	float feedbackI = state[2];
  	float feedbackQ = state[3];
  	ncoOut[0] = state[4];
  	float trigOffset = state[5];

		float errorI, errorQ, errorD, trigArg;

    for (int k = 0; k < pllIn.size(); k++){
      errorI = pllIn[k] * (+feedbackI);  //complex conjugate of the
  		errorQ = pllIn[k] * (-feedbackQ);  //feedback complex exponential

  		//four-quadrant arctangent discriminator for phase error detection
  		errorD = std::atan2(errorQ, errorI);

  		//loop filter
  		integrator = integrator + (Ki*errorD);

  		//update phase estimate
  		phaseEst = phaseEst + (Kp*errorD) + integrator;

  		//internal oscillator
  		trigOffset += 1.0;
  		trigArg = (2*PI*(freq/Fs)*(trigOffset) + phaseEst);
  		feedbackI = std::cos(trigArg);
  		feedbackQ = std::sin(trigArg);
  		ncoOut[k+1] = std::cos(trigArg*ncoScale + phaseAdjust);
    }
    state[0] = integrator;
		state[1] = phaseEst;
		state[2] = feedbackI;
		state[3] = feedbackQ;
		state[4] = ncoOut[-1];
		state[5] = trigOffset;
}

void mixer(std::vector<float> &recoveredStereo, std::vector<float> &channel_filt, std::vector<float> &mixedAudio)
{
  mixedAudio.clear();
  mixedAudio.resize(channel_filt.size());

  for (int i = 0; i < channel_filt.size(); i++){
    mixedAudio[i] = 2 * recoveredStereo[i] * channel_filt[i]	//this would have both the +ve and -ve part of the cos combined, we need to keep the -ve part and filter it
  }
}
<<<<<<< HEAD

std::vector<float> stereo_path_main(int mode){
	// TIMING VARIABLES
  auto start_overall_time = std::chrono::high_resolution_clock::now();
  auto start_time = std::chrono::high_resolution_clock::now();
  auto stop_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> STEREO_run_time;
  std::chrono::duration<double, std::milli> RF_run_time;

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
	float audio_Fs = 48000;
	int audio_taps = 151;

	std::vector<float> rf_coeff, carrier_coeff, channel_coeff, audio_coeff;

	low_pass_coeff(rf_Fs, rf_Fc, rf_taps, rf_coeff);
	std::cerr << "test 1" << "\n";
	band_pass_coeff(18500, 19500, rf_Fs/rf_decim, audio_taps, carrier_coeff);
	std::cerr << "test 2" << "\n";
	band_pass_coeff(22000, 54000, rf_Fs/rf_decim, audio_taps, channel_coeff);
	std::cerr << "test 3" << "\n";
	low_pass_coeff((rf_Fs/rf_decim)*audio_exp, audio_Fc, audio_taps*audio_exp, audio_coeff);
	std::cerr << "test 4" << "\n";

	float block_size = 1024*rf_decim*audio_decim*2;

  // filtered data
  std::vector<float> i_filt;
	std::vector<float> q_filt;

  // demodulation variables
	std::vector<float> demod_state;
	std::vector<float> IQ_demod;
	demod_state.resize(2, 0.0);

	std::vector<float> filt_block, carrier_block, channel_block, mono_block, pll_block;
	filt_block.resize(carrier_coeff.size());		//For Stereo Processing
	carrier_block.resize(carrier_coeff.size());	//For Stereo Carrier Recovery
	channel_block.resize(channel_coeff.size());	//For Stereo Channel Extraction
	mono_block.resize(audio_coeff.size());		//For Monopath

	std::vector<float> left_block, right_block;
	left_block.resize(1024);
	right_block.resize(1024);

	pll_block.resize(6);

	pll_block[0] = 0.0;
	pll_block[1] = 0.0;
	pll_block[2] = 1.0;
	pll_block[3] = 0.0;
	pll_block[4] = 1.0;
	pll_block[5] = 0.0;

	// state saving variables for I and Q samples convolution
	std::vector<float> state_i_lpf_100k;
	std::vector<float> state_q_lpf_100k;
	state_i_lpf_100k.resize(rf_coeff.size() - 1, 0.0);
	state_q_lpf_100k.resize(rf_coeff.size() - 1, 0.0);

	//audio buffer that stores all the audio blocks
	std::vector<float> audio_data, stereo_data, left_data, right_data;
	std::vector<float> carrier_filt, dummy_fm, channel_filt, down_carrier, down_channel, recoveredStereo, down_stereo;

	// STEP 1: IQ samples demodulation
	std::vector<float> i_data(block_size / 2);
	std::vector<float> q_data(block_size / 2);

	carrier_filt.resize(i_data.size());
	channel_filt.resize(i_data.size());

	std::vector<std::vector<float>> complete_data;

  // decipher each block
	for(unsigned int block_id = 0; ; block_id++) {
		std::vector<float> block_data(block_size);
    readStdinBlockData(block_size, block_data);
    if ((std::cin.rdstate()) != 0) {
      std::cerr << "End of input stream reached" << std::endl;

      // timing analysis
      stop_time = std::chrono::high_resolution_clock::now();
		  std::chrono::duration<double, std::milli> OVERALL_run_time = stop_time-start_overall_time;
		  std::cerr << "OVERALL RUNTIME: " << OVERALL_run_time.count() << " ms" << "\n";
      std::cerr << "RF DOWNSAMPLE RUNTIME: " << RF_run_time.count() << " ms" << "\n";
      std::cerr << "TOTAL MONOPATH RUNTIME: " << STEREO_run_time.count() << " ms" << "\n";
      exit(1);
    }

		for (int k = 0; k < block_size / 2; k++) {
      i_data[k] = block_data[2*k];
      q_data[k] = block_data[2*k + 1];
    }

    // filter out IQ data with convolution
    start_time = std::chrono::high_resolution_clock::now();
    // take downsampled filtered IQ data
		std::vector<float>i_ds;
    i_ds.resize(i_filt.size()/rf_decim);
		std::vector<float>q_ds;
    q_ds.resize(q_filt.size()/rf_decim);

		ds_block_conv(i_filt,i_data,rf_coeff,state_i_lpf_100k,rf_decim,i_ds);
		ds_block_conv(q_filt,q_data,rf_coeff,state_q_lpf_100k,rf_decim,q_ds);

    // timing analysis
    stop_time = std::chrono::high_resolution_clock::now();
    RF_run_time += stop_time-start_time;

		// perform demodulation on IQ data
		IQ_demod = fmDemod(i_ds, q_ds, demod_state);
    //std::cerr << "test1\n";
		// STEP 2: Mono path
    start_time = std::chrono::high_resolution_clock::now();
    //std::vector<float> audio_block = mono_path(mode, IQ_demod, audio_coeff, audio_state, audio_decim, audio_exp);

		//Stereo Carrier Recovery: Bandpass -> PLL -> Numerically Controlled
		const std::vector<float> x;
		std::cerr << "test 8" << "\n";
		ds_block_conv(carrier_filt, IQ_demod, carrier_coeff, carrier_block, rf_decim, down_carrier);
		std::cerr << "test 9" << "\n";
		fmPll(carrier_filt, recoveredStereo, pll_block, 19e3, rf_Fs/rf_decim, 2.0, 0.0, 0.01);
		std::cerr << "test 5" << "\n";

    //Stereo Channel Extraction: Bandpass
		ds_block_conv(channel_filt, IQ_demod, channel_coeff, channel_block, rf_decim, down_channel);

    //Stereo Processing: Mixer -> Digital filtering (Lowpass -> down sample) -> Stereo Combiner
		std::vector<float> mixedAudio, stereo_filt, stereo_block;
		mixer(recoveredStereo, channel_filt, mixedAudio);
		std::cerr << "test 6" << "\n";

		ds_block_conv(stereo_filt, mixedAudio, audio_coeff, filt_block, rf_decim, down_stereo);

		for (int i = 0; i < stereo_filt.size(); i+=audio_decim){
			stereo_block.push_back(stereo_filt[i]);
		}

		//extract the mono audio data
		std::vector<float> audio_filt, audio_block;

		// state saving variable for audio data convolution
		std::vector<float> audio_state;
		audio_state.resize(audio_coeff.size() - 1, 0.0);

    ds_block_conv(audio_filt, IQ_demod, audio_coeff, audio_state, audio_decim, audio_block);
		audio_block.reserve(audio_filt.size()/audio_decim);

		std::vector<float> index;

		for (int i = 0; i < audio_block.size(); i++){
			index.clear();
			left_block[i] = audio_block[i] + stereo_block[i];
			right_block[i] = audio_block[i] - stereo_block[i];
			index.push_back(left_block[i]);
			index.push_back(right_block[i]);
			complete_data.push_back(index);
		}

		//left_data.insert(left_data.end(), left_block.begin(), left_block.end());
		//right_data.insert(right_data.end(), right_data.begin(), right_data.end());

    // timing analysis
    stop_time = std::chrono::high_resolution_clock::now();
    STEREO_run_time += stop_time-start_time;

    // STEP 3: prepare audio data for output
    std::vector<std::vector<short int>> audio_data;
		std::vector<short int> complete_block;

    for (unsigned int k = 0; k < complete_data.size(); k++) {
      if(std::isnan(complete_data[k][0])){
				std::vector<short int> temp;
				temp.push_back(0);
				temp.push_back(0);
        audio_data.push_back(temp);
      }
      else{
				complete_block.clear();
				complete_block.push_back(static_cast<short int>(complete_data[k][0] * 16384));
				complete_block.push_back(static_cast<short int>(complete_data[k][1] * 16384));
				audio_data.push_back(complete_block);
      }
    }
    fwrite(&audio_data[0], sizeof(short int), audio_data.size(), stdout);

	}
	std::cerr << "test 7" << "\n";

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
  stereo_path_main(mode);

	return 0;
}
=======
>>>>>>> 053f26627ab76c199cdca6eb75671c38cbf867f3
