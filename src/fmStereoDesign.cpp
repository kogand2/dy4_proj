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

//Low pass filter coefficients
void low_pass_coeff(float Fs, float Fc, int num_taps, std::vector<float> &h)
{
	// allocate memory for the impulse response
	h.clear();
	h.resize(num_taps, 0.0);

	float norm_co = Fc/(Fs/2);

	// the rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
	for (int i = 0; i < num_taps; i++){
			if (i == (num_taps-1)/2)
		    h[i] = norm_co;
			else
			  h[i] = norm_co * (std::sin(PI*norm_co*(i-((num_taps-1)/2))))/(PI*norm_co*(i-((num_taps-1)/2)));
			h[i] = h[i] * pow(std::sin(i*PI/num_taps), 2);
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

// block convolution function (with downsampling)
void ds_block_conv(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &state, int rf_decim, std::vector<float> &down)
{
  	// allocate memory for the output (filtered) data
  	y.clear();
  	y.resize(x.size(), 0.0); // y of size i_data

    // clear downsampled output
    down.clear();

    // only compute the values we need (because of downsampling)
  	for (int n = 0; n < y.size(); n += rf_decim){
  		for (int k = 0; k < h.size(); k++){
  			if (n-k >= 0){
  				y[n] += h[k] * x[n - k];
        }
        else{
  				y[n] += h[k] * state[state.size() + n - k];
        }
      }
      down.push_back(y[n]);
    }

    int index = x.size() - h.size() + 1;
  	state = std::vector<float>(x.begin() + index, x.end());
  }

void fmPll(std::vector<float> &pllIn, std::vector<float> &ncoOut, std::vector<float> &new_state, float freq, float Fs, float ncoScale = 2.0, float phaseAdjust = 0.0, float normBandwidth = 0.01, std::vector<float> &state)
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
    new_state = [integrator, phaseEst, feedbackI, feedbackQ, ncoOut[-1], trigOffset];
}

void mixer(std::vector<float> &recoveredStereo, std::vector<float> &channel_filt, std::vector<float> &mixedAudio)
{
  mixedAudio.clear();
  mixedAudio.resize(channel_filt.size());

  for (int i = 0; i < channel_filt.size(); i++){
    mixedAudio[i] = 2 * recoveredStereo[i] * channel_filt[i];	//this would have both the +ve and -ve part of the cos combined, we need to keep the -ve part and filter it
  }
}

void readStdinBlockData(unsigned int num_samples, std::vector<float> &block_data){
  std::vector<char> raw_data(num_samples);
  std::cin.read(reinterpret_cast<char*>(&raw_data[0]), num_samples*sizeof(char));
  for(int k = 0; k < (int)num_samples; k++) {
    block_data[k] = float(((unsigned char)raw_data[k] - 128)/ 128.0);
  }
}


std::vector<float> stereo_path_main(){
	// TIMING VARIABLES
  auto start_overall_time = std::chrono::high_resolution_clock::now();
  auto start_time = std::chrono::high_resolution_clock::now();
  auto stop_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> MONO_run_time;
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
	band_pass_coeff(18500, 19500, rf_Fs/rf_decim, audio_taps, carrier_coeff);
	band_pass_coeff(22000, 54000, rf_Fs/rf_decim, audio_taps, channel_coeff);
	low_pass_coeff((rf_Fs/rf_decim)*audio_exp, audio_Fc, audio_taps*audio_exp, audio_coeff);

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

	std::vector<float> prevStereo, prevCarrier;
	prevStereo.resize(5121);
	prevCarrier.resize(5120);

	// state saving variables for I and Q samples convolution
	std::vector<float> state_i_lpf_100k;
	std::vector<float> state_q_lpf_100k;
	state_i_lpf_100k.resize(rf_coeff.size() - 1, 0.0);
	state_q_lpf_100k.resize(rf_coeff.size() - 1, 0.0);

	//audio buffer that stores all the audio blocks
	std::vector<float> audio_data, stereo_data, left_data, right_data;

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
      std::cerr << "TOTAL MONOPATH RUNTIME: " << MONO_run_time.count() << " ms" << "\n";
      exit(1);
    }

		// STEP 1: IQ samples demodulation
		std::vector<float> i_data(block_size / 2);
    std::vector<float> q_data(block_size / 2);

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
		std::vector<float> carrier_filt, down_carrier, down_channel, recoveredStereo, new_state, down_stereo;
		const std::vector<float> x;

		ds_block_conv(carrier_filt, dummy_fm, carrier_coeff, carrier_block, rf_decim, down_carrier);
		fmPll(carrier_filt, recoveredStereo, new_state, 19e3, rf_Fs/rf_decim, 2.0, 0.0, 0.01, pll_block);

    //Stereo Channel Extraction: Bandpass
		ds_block_conv(channel_filt, dummy_fm, channel_coeff, channel_block, rf_decim, down_channel);

    //Stereo Processing: Mixer -> Digital filtering (Lowpass -> down sample) -> Stereo Combiner
		std::vector<float> mixedAudio, stereo_filt, stereo_block;
		mixer(recoveredStereo, channel_filt, mixedAudio);

		ds_block_conv(stereo_filt, mixedAudio, audio_coeff, filt_block, rf_decim, down_stereo);

		for (int i = 0; i < stereo_filt.size(); i+=audio_decim){
			stereo_block.push_back(stereo_filt[i)];
		}

		//extract the mono audio data
		std::vector<float> audio_filt, audio_block;

		// state saving variable for audio data convolution
		std::vector<float> audio_state;
		audio_state.resize(audio_coeff.size() - 1, 0.0);

    ds_block_conv(audio_filt, IQ_demod, audio_coeff, audio_state, audio_decim, audio_block);
		audio_block.reserve(audio_filt.size()/audio_decim);

		for (int i = 0; i < audio_block.size(); i++){
			left_block[i] = audio_block[i] + stereo_block[i];
			right_block[i] = audio_block[i] - stereo_block[i];
		}

		left_data.insert(left_data.end(), left_block.begin(), left_block.end());
		right_data.insert(right_data.end(), right_data.begin(), right_data.end());

    // timing analysis
    stop_time = std::chrono::high_resolution_clock::now();
    STEREO_run_time += stop_time-start_time;

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
