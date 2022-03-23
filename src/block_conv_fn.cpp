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

// Filter coefficient function
void low_pass_coeff(float Fs, float Fc, int num_taps, std::vector<float> &h)
{
	// allocate memory for the impulse response
	//h.clear();
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

// convolution with no downsampling
void state_block_conv(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &state)
{
	// allocate memory for the output (filtered) data
	y.clear();
	y.resize(x.size(), 0.0);

	// implement block processing algorithm discussed in lecture and used in python
	for (int n = 0; n < y.size(); n++){
		for (int k = 0; k < h.size(); k++){
			if (n-k >= 0){
				y[n] += h[k] * x[n - k];
        //std::cerr << "y[" << n << "] += h[" << k << "] * x[" << n - k << "]\n";
			}else{
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
	// allocate memory for the output (filtered) data
  auto start_time = std::chrono::high_resolution_clock::now();
	//y.clear();
	//y.resize(x.size(), 0.0); // y of size i_data
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
	// allocate memory for the output (filtered) data
  auto start_time = std::chrono::high_resolution_clock::now();
	//y.clear();
	//y.resize(x.size()*audio_exp, 0.0); // y of size i_data
  int count = 0;
  int phase, x_index;

  auto stop_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> SFT_run_time = stop_time-start_time;
  std::cerr << "PREP RUNTIME: " << SFT_run_time.count() << " ms" << "\n";

  start_time = std::chrono::high_resolution_clock::now();
  // only compute the values we need (because of downsampling)
	for (int n = 0; n < x.size()*audio_exp; n += audio_decim){
    phase = n % audio_exp;
    x_index = (n - phase) / audio_exp;
    y[n] = 0;
		for (int k = phase; k < h.size(); k += audio_exp){
			if (x_index >= 0){
				y[n] += audio_exp * h[k] * x[x_index];
        //std::cerr << "y[" << n << "] += h[" << k << "] * x[" << x_index << "]\n";
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
  std::cerr << "FOR LOOP RUNTIME: " << DFT_run_time.count() << " ms" << "\n";

  start_time = std::chrono::high_resolution_clock::now();

  std::cerr << "rs: " << count << "\n";
  int index = x.size() - h.size()/audio_exp + 1;
	state = std::vector<float>(x.begin() + index, x.end());

  stop_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> NFT_run_time = stop_time-start_time;
  std::cerr << "STATE SAVING RUNTIME: " << NFT_run_time.count() << " ms" << "\n";
}
