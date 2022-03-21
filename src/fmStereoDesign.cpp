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
    h[i] = 2.5 * h[i] * pow(std::sin(i*PI/num_taps), 2);
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
    mixedAudio[i] = 2 * recoveredStereo[i] * channel_filt[i]	//this would have both the +ve and -ve part of the cos combined, we need to keep the -ve part and filter it
  }
}
