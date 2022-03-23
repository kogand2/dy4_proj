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

void readStdinBlockData(unsigned int num_samples, std::vector<float> &block_data){
  std::vector<char> raw_data(num_samples);
  std::cin.read(reinterpret_cast<char*>(&raw_data[0]), num_samples*sizeof(char));
  for(int k = 0; k < (int)num_samples; k++) {
    block_data[k] = float(((unsigned char)raw_data[k] - 128)/ 128.0);
  }
}

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
    h[i] = h[i] * pow(std::sin(i*PI/num_taps), 2);
  }
}

void all_pass_coeff(std::vector<float> &y,std::vector<float> &x, std::vector<float> &all_pass_state)
{
  y.clear();
  y = all_pass_state;
  std::vector<float> curr_state = std::vector<float>(x.begin(), x.end() - all_pass_state.size());
  y.insert(y.end(), curr_state.begin(), curr_state.end());

  all_pass_state = std::vector<float>(x.begin() + x.size() - all_pass_state.size(), x.end());
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
		state[4] = ncoOut[ncoOut.size() - 1];
		state[5] = trigOffset;
}

void mixer(std::vector<float> &recoveredStereo, std::vector<float> &channel_filt, std::vector<float> &mixedAudio)
{
  mixedAudio.clear();
//  mixedAudio.resize(channel_filt.size());

  for (int i = 0; i < channel_filt.size(); i++){
    mixedAudio.push_back(2 * recoveredStereo[i] * channel_filt[i]);	//this would have both the +ve and -ve part of the cos combined, we need to keep the -ve part and filter it
  }
  std::cerr << mixedAudio.size() << "\n";
  std::cerr << channel_filt.size() << "\n";
  //exit(0);
}

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
	int audio_taps = 101;

	std::vector<float> rf_coeff, carrier_coeff, channel_coeff, audio_coeff;

  // filter coefficients for RF, Stereo, Mono
	low_pass_coeff(rf_Fs, rf_Fc, rf_taps, rf_coeff);
	band_pass_coeff(18500, 19500, rf_Fs/rf_decim, audio_taps, carrier_coeff);
	band_pass_coeff(22000, 54000, rf_Fs/rf_decim, audio_taps, channel_coeff);
	low_pass_coeff((rf_Fs/rf_decim)*audio_exp, audio_Fc, audio_taps*audio_exp, audio_coeff);

	float block_size = 102400;

  // demodulation variables
	std::vector<float> demod_state;
	std::vector<float> IQ_demod;
	demod_state.resize(2, 0.0);

	std::vector<float> stereo_state, carrier_state, channel_state, mono_state, pll_state;
	stereo_state.resize(carrier_coeff.size() - 1);		//For Stereo Processing
	carrier_state.resize(carrier_coeff.size() - 1);	//For Stereo Carrier Recovery
	channel_state.resize(channel_coeff.size() - 1);	//For Stereo Channel Extraction
	mono_state.resize(audio_coeff.size() - 1);		  //For Monopath

  // might have to change block size
	std::vector<float> left_block, right_block, delay_block;
  std::vector<float> mono_input;
	left_block.resize(1024);
	right_block.resize(1024);
  delay_block.resize((audio_taps-1)/2);

	pll_state.resize(6);

	pll_state[0] = 0.0;
	pll_state[1] = 0.0;
	pll_state[2] = 1.0;
	pll_state[3] = 0.0;
	pll_state[4] = 1.0;
	pll_state[5] = 0.0;

	// state saving variables for I and Q samples convolution
	std::vector<float> state_i_lpf_100k;
	std::vector<float> state_q_lpf_100k;
	state_i_lpf_100k.resize(rf_coeff.size() - 1, 0.0);
	state_q_lpf_100k.resize(rf_coeff.size() - 1, 0.0);

	//audio buffer that stores all the audio blocks
	std::vector<float> audio_data, stereo_data, left_data, right_data;
	std::vector<float> carrier_filt, dummy_fm, channel_filt, down_carrier, down_channel, recoveredStereo;

	// STEP 1: IQ samples demodulation
	std::vector<float> i_data(block_size / 2);
	std::vector<float> q_data(block_size / 2);

  // filtered IQ data
  std::vector<float> i_filt;
  i_filt.resize(i_data.size());
  std::vector<float> q_filt;
  q_filt.resize(q_data.size());

  // take downsampled filtered IQ data
  std::vector<float>i_ds;
  i_ds.resize(i_filt.size()/rf_decim);
  std::vector<float>q_ds;
  q_ds.resize(q_filt.size()/rf_decim);

	carrier_filt.resize(i_data.size());
	channel_filt.resize(i_data.size());

  // filtered audio data
  std::vector<float> audio_filt;
  audio_filt.resize(audio_exp*(i_ds.size() - 1));

  // downsampled filtered audio data
  std::vector<float> audio_block;
  audio_block.resize(audio_filt.size()/audio_decim);

  // state saving variable for audio data convolution
  std::vector<float> audio_state;
  audio_state.resize(audio_coeff.size() - 1);

  std::vector<float> mixedAudio, stereo_filt, stereo_block;

  mixedAudio.resize(channel_filt.size());
  stereo_filt.resize((mixedAudio.size()-1)*audio_exp);


	std::vector<std::vector<float>> complete_data;
	//std::cerr << "test 5" << "\n";
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

    std::cerr << "Processing Block " << block_id << std::endl;

		for (int k = 0; k < block_size / 2; k++) {
      i_data[k] = block_data[2*k];
      q_data[k] = block_data[2*k + 1];
    }

    // filter out IQ data with convolution
    start_time = std::chrono::high_resolution_clock::now();
    i_ds.clear();
    q_ds.clear();
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
		//std::cerr << "test 8" << "\n";
    down_carrier.clear();
		state_block_conv(carrier_filt, IQ_demod, carrier_coeff, carrier_state);

		fmPll(carrier_filt, recoveredStereo, pll_state, 19e3, rf_Fs/rf_decim, 2.0, 0.0, 0.01);

    //Stereo Channel Extraction: Bandpass
    down_channel.clear();
		state_block_conv(channel_filt, IQ_demod, channel_coeff, channel_state);

    //Stereo Processing: Mixer -> Digital filtering (Lowpass -> down sample) -> Stereo Combiner
		mixer(recoveredStereo, channel_filt, mixedAudio);



    stereo_block.clear();

    if (mode == 0 || mode == 1)
		  ds_block_conv(stereo_filt, mixedAudio, audio_coeff, stereo_state, audio_decim, stereo_block);
    else
      rs_block_conv(stereo_filt, mixedAudio, audio_coeff, stereo_state, audio_decim, audio_exp, stereo_block);

    all_pass_coeff(mono_input, IQ_demod, delay_block);


    audio_block.clear();
    if (mode == 0 || mode == 1)
      ds_block_conv(audio_filt, mono_input, audio_coeff, audio_state, audio_decim, audio_block);
    else
      rs_block_conv(audio_filt, mono_input, audio_coeff, audio_state, audio_decim, audio_exp, audio_block);

    //printRealVector(audio_block);
    //exit(0);

		std::vector<float> index;

    complete_data.clear();
		for (int i = 0; i < audio_block.size(); i++){
      //std::cerr << "test 9" << "\n";
			float left_block = audio_block[i] + stereo_block[i];
			float right_block = audio_block[i] - stereo_block[i];

      //std::cerr << "Left Block [ " << i << " ] = " <<  left_block << "\n";
      //std::cerr << "Right Block [ " << i << " ] = " <<  right_block << "\n";

			index.push_back(left_block);
			index.push_back(right_block);
      //std::cerr << "test 11" << "\n";
			complete_data.push_back(index);
      index.clear();
		}
    //std::cerr << "test2\n";

    //std::cerr << "test3\n";
    //std::cerr << "test 9" << "\n";
		//left_data.insert(left_data.end(), left_block.begin(), left_block.end());
		//right_data.insert(right_data.end(), right_data.begin(), right_data.end());

    // timing analysis
    stop_time = std::chrono::high_resolution_clock::now();
    STEREO_run_time += stop_time-start_time;

    // STEP 3: prepare audio data for output
    std::vector<short int> audio_data;
		std::vector<short int> complete_block;

    //std::cerr << "test 9\n";
    for (unsigned int k = 0; k < complete_data.size(); k++) {
      if(std::isnan(complete_data[k][0])){
        audio_data.push_back(0);
        audio_data.push_back(0);
      }
      else{
				//complete_block.clear();
				audio_data.push_back(static_cast<short int>(complete_data[k][0] * 16384));
				audio_data.push_back(static_cast<short int>(complete_data[k][1] * 16384));
        //std::cerr << "Complete Block [ " << k << "] " << complete_data[k][0] << ", " << complete_data[k][1] << "\n";

        //audio_data.push_back(complete_block);
      }
    }
    /*if (block_id == 2)
    {
      exit(0);
    }*/
    //break;
    //std::cerr << "test 10" << "\n";
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
  stereo_path_main(mode);

	return 0;
}
