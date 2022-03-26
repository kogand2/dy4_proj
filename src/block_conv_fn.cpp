#include "block_conv_fn.h"
#include "dy4.h"
#include "iofunc.h"

// testing time complexity
#include <chrono>

// Custom Demodulation Function
std::vector<float> fmDemod(std::vector<float> I, std::vector<float> Q, std::vector<float> &dummy_state) {
  std::vector<float> fm_demod;
  std::vector<float> I_new;
  std::vector<float> Q_new;
  float I_prev;
  float Q_prev;
  float I_der;
  float Q_der;
  int j;

  fm_demod.resize(I.size() - 1, 0.0);

  I_prev = dummy_state[0];
  Q_prev = dummy_state[1];

  I_new.push_back(I_prev);
  Q_new.push_back(Q_prev);

  I_new.insert( I_new.end(), I.begin(), I.end() );
  Q_new.insert( Q_new.end(), Q.begin(), Q.end() );

  for (j = 0; j < I.size() - 1; j++){
    Q_der = Q_new[j+1] - Q_new[j];
    I_der = I_new[j+1] - I_new[j];

    fm_demod[j] = Q_der*I_new[j+1] - I_der*Q_new[j+1];

    if (pow(Q_new[j+1], 2.0) + pow(I_new[j+1], 2.0) != 0){
      fm_demod[j] = fm_demod[j]/(pow(Q_new[j+1], 2.0) + pow(I_new[j+1], 2.0));
    }else{
      fm_demod[j] = 0;
    }
  }

  dummy_state.clear();

  dummy_state[0] = I[j];
  dummy_state[1] = Q[j];

  return fm_demod;
}

// LPF coefficient function
void low_pass_coeff(float Fs, float Fc, int num_taps, std::vector<float> &h)
{
	// allocate memory for the impulse response
	h.resize(num_taps, 0.0);

	float norm_co = Fc/(Fs/2);

	for (int i = 0; i < num_taps; i++){
			if (i == (num_taps-1)/2)
		    h[i] = norm_co;
			else
			  h[i] = norm_co * (std::sin(PI*norm_co*(i-((num_taps-1)/2))))/(PI*norm_co*(i-((num_taps-1)/2)));

      h[i] = h[i] * pow(std::sin(i*PI/num_taps), 2);
	}
}

// BPF coefficient function
void band_pass_coeff(float Fb, float Fe, float Fs, int num_taps, std::vector<float> &h)
{
  h.clear();
  h.resize(num_taps, 0.0);

  float normCenter = ((Fe + Fb)/2)/(Fs/2);
	float normPass = (Fe - Fb)/(Fs/2);

  for (int i = 0; i < num_taps; i++){
    if (i == (num_taps-1)/2)
      h[i] = normPass;
    else
      h[i] = normPass * (std::sin(PI*(normPass/2)*(i-((num_taps-1)/2))))/(PI*(normPass/2)*(i-((num_taps-1)/2)));

    h[i] = h[i] * std::cos(i*PI*normCenter);
    h[i] = h[i] * pow(std::sin(i*PI/num_taps), 2);
  }
}

// APF coefficient function
void all_pass_coeff(std::vector<float> &y,std::vector<float> &x, std::vector<float> &all_pass_state)
{
  y.clear();
  y = all_pass_state;
  std::vector<float> curr_state = std::vector<float>(x.begin(), x.end() - all_pass_state.size());
  y.insert(y.end(), curr_state.begin(), curr_state.end());

  all_pass_state = std::vector<float>(x.begin() + x.size() - all_pass_state.size(), x.end());
}

// FM PLL function
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

// mixer function
void mixer(std::vector<float> &recoveredStereo, std::vector<float> &channel_filt, std::vector<float> &mixedAudio)
{
  mixedAudio.clear();
//  mixedAudio.resize(channel_filt.size());

  for (int i = 0; i < channel_filt.size(); i++){
    mixedAudio.push_back(2 * recoveredStereo[i] * channel_filt[i]);	//this would have both the +ve and -ve part of the cos combined, we need to keep the -ve part and filter it
  }
}

std::vector<float> fmDemodArctan(std::vector<float> I, std::vector<float> Q, float &prev_phase)
{
  std::vector<float> fm_demod;
  fm_demod.resize(I.size());
  float current_phase;
  float phase_diff;
  for (int k = 0; k < I.size(); k++){
    current_phase = std::atan2(Q[k], I[k]);

    phase_diff = current_phase - prev_phase;

    if (fabs(phase_diff) != PI){
      phase_diff = fmod(phase_diff + PI, 2*PI);


      if (phase_diff < 0){
        phase_diff += 2*PI;
      }

      phase_diff -= PI;
    }

    fm_demod[k] = prev_phase + phase_diff;

    prev_phase = current_phase;

  }
  return fm_demod;

}

// convolution with no downsampling
void state_block_conv(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &state)
{
  // FOR NOW
  y.clear();
	y.resize(x.size(), 0.0);

	// implement block processing algorithm discussed in lecture and used in python
	for (int n = 0; n < y.size(); n++){
    y[n] = 0;
		for (int k = 0; k < h.size(); k++){
			if (n-k >= 0){
				y[n] += h[k] * x[n - k];
        //std::cerr << "y[" << n << "] += h[" << k << "] * x[" << n - k << "]\n";
			}
      else{
				y[n] += h[k] * state[state.size() + n - k];
        //std::cerr << "y[" << n << "] += h[" << k << "] * state[" << state.size() +  n - k << "]\n";
      }
    }
  }

  int index = x.size() - h.size() + 1;
	state = std::vector<float>(x.begin() + index, x.end());

}

// block convolution function (with downsampling)
void ds_block_conv(std::vector<float> &y, const std::vector<float> x, const std::vector<float> h, std::vector<float> &state, int rf_decim, std::vector<float> &down)
{
  auto start_time = std::chrono::high_resolution_clock::now();
  auto stop_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> SFT_run_time = stop_time-start_time;
  //std::cerr << "PREP RUNTIME: " << SFT_run_time.count() << " ms" << "\n";

  int count = 0;
  //start_time = std::chrono::high_resolution_clock::now();
  // only compute the values we need (because of downsampling)
	for (int n = 0; n < x.size(); n += rf_decim){
    y[n] = 0;
		for (int k = 0; k < h.size(); k++){
			if (n-k >= 0){
				y[n] += h[k] * x[n - k];
      }
      else{
				y[n] += h[k] * state[state.size() + n - k];
      }
      count++;
    }
    down.push_back(y[n]);
  }

  stop_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> DFT_run_time = stop_time-start_time;
  //std::cerr << "FOR LOOP RUNTIME: " << DFT_run_time.count() << " ms" << "\n";

  start_time = std::chrono::high_resolution_clock::now();

  //std::cerr << "ds: " << count << "\n";
  int index = x.size() - h.size() + 1;
	state = std::vector<float>(x.begin() + index, x.end());

  stop_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> NFT_run_time = stop_time-start_time;
  //std::cerr << "STATE SAVING RUNTIME: " << NFT_run_time.count() << " ms" << "\n";
}

// block convolution function (with resampling)
void rs_block_conv(std::vector<float> &y, const std::vector<float> x, const std::vector<float> h, std::vector<float> &state, int audio_decim, int audio_exp, std::vector<float> &down)
{
  auto start_time = std::chrono::high_resolution_clock::now();
  int count = 0;
  int phase, x_index;

  auto stop_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> SFT_run_time = stop_time-start_time;
  //std::cerr << "PREP RUNTIME: " << SFT_run_time.count() << " ms" << "\n";

  start_time = std::chrono::high_resolution_clock::now();
  // only compute the values we need (because of downsampling)
	for (int n = 0; n < x.size()*audio_exp; n += audio_decim){

    phase = n % audio_exp;

    x_index = (n - phase) / audio_exp;

    y[n] = 0;

		for (int k = phase; k < h.size(); k += audio_exp){
			if (x_index >= 0){

        //std::cerr << "y[" << n << "] += h[" << k << "] * x[" << x_index << "]\n";
				y[n] += audio_exp * h[k] * x[x_index];

      }
      else{

				y[n] += audio_exp * h[k] * state[state.size() + x_index];
        //std::cerr << "y[" << n << "] += h[" << k << "] * state[" << state.size() +  x_index << "]\n";
      }
      x_index--;
      count++;
    }
    down.push_back(y[n]);

  }


  stop_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> DFT_run_time = stop_time-start_time;
  //std::cerr << "FOR LOOP RUNTIME: " << DFT_run_time.count() << " ms" << "\n";

  start_time = std::chrono::high_resolution_clock::now();

  //std::cerr << "rs: " << count << "\n";
  int index = x.size() - h.size()/audio_exp + 1;
	state = std::vector<float>(x.begin() + index, x.end());

  stop_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> NFT_run_time = stop_time-start_time;
  //std::cerr << "STATE SAVING RUNTIME: " << NFT_run_time.count() << " ms" << "\n";
}
