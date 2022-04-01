#include "block_conv_fn.h"
#include "dy4.h"
#include "iofunc.h"
#include <math.h>

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

float unwrap(float previous_angle, float new_angle) {
    float d = new_angle - previous_angle;
    d = d > M_PI ? d - 2 * M_PI : (d < -M_PI ? d + 2 * M_PI : d);
    return previous_angle + d;
}

std::vector<float> fmDemodArctan(std::vector<float> I, std::vector<float> Q, float &prev_phase)
{
  std::vector<float> fm_demod;
  fm_demod.resize(I.size());
  float current_phase;
  for (int k = 0; k < I.size(); k++){
    current_phase = std::atan2(Q[k], I[k]);

    fm_demod[k] = unwrap(prev_phase, current_phase) - prev_phase;

    prev_phase = current_phase;

  }
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

// RRC coefficient function
void rrc_coeff(float Fs, int num_taps, std::vector<float> &h)
{
  h.clear();
  h.resize(num_taps, 0.0);

  float T_symbol = 1/2375.0;
  float beta = 0.90;
  float t;

  for (int i = 0; i < num_taps; i++){
    t = (i - num_taps/2.0)/Fs;
    if (t == 0.0){
      h[i] = 1.0 + beta*((4/PI) - 1);
    }else if (t == -T_symbol/(4*beta) || t == T_symbol/(4*beta)){
      h[i] = (beta/sqrt(2))*(((1 + 2/PI)*(sin(PI/(4*beta)))) + ((1 - 2/PI)*(cos(PI/(4*beta)))));
    }else{
      h[i] = (sin(PI*t*(1 - beta)/T_symbol) + 4*beta*(t/T_symbol)*cos(PI*t*(1+beta)/T_symbol))/(PI*t*(1 - (4*beta*t/T_symbol)*(4*beta*t/T_symbol))/T_symbol);
    }
  }
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

  for (int i = 0; i < channel_filt.size(); i++){
    mixedAudio.push_back(2 * recoveredStereo[i] * channel_filt[i]);	//this would have both the +ve and -ve part of the cos combined, we need to keep the -ve part and filter it
  }
}

// square non-linearity function
std::vector<float> sq_non_linearity(std::vector<float> signalIn)
{
  std::vector<float> carrier_input;
  for (int i = 0; i < signalIn.size(); i++)
    carrier_input.push_back(signalIn[i] * signalIn[i]);

  return carrier_input;
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

void CDR(std::vector<float> signalIn, int interval, int &initial, std::vector<int> &sample_idxs, std::vector<int> &man_encoding){

  sample_idxs.clear();	// x vals (just for testing)
  man_encoding.clear();			// y vals (manchester_values)

  //# processing first block, find best initial point to start sampling
  if (initial == -1){
    float max = 0.0;
    float curr = 0.0;
    int best_init = 0;

    for (int init = 0; init < interval; init++){
      // reset for next vals
      curr = 0.0;
      for (int i = init; i < signalIn.size(); i += interval){
        curr += abs(signalIn[i]);
      }
      if (max < curr){
        max = curr;
        initial = init;
      }
    }
  }

	// processing rest of blocks, finding values of samples
	for (int i = initial; i < signalIn.size(); i += interval){
		sample_idxs.push_back(i);
		if (signalIn[i] > 0){	// get a HI
			man_encoding.push_back(1);
		}else{					// get a LO
			man_encoding.push_back(0);
    }
  }
}

int diff_decoding(std::vector<int> manchester_values, std::vector<int> cdr_state, std::vector<int> &decoded_bits, int &initial){
  std::vector<int> bits;
  decoded_bits.clear();
  bool is_first_block = false;

  if (initial == -1){
    is_first_block = true;
    int ind1 = 0;
    int ind2 = 0;

    for (int i = 0; i < manchester_values.size() - 1; i += 2){
      if (manchester_values[i] == 1 && manchester_values[i+1] == 1){ // HH (ignore)
        ind1++;
      }else if (manchester_values[i] == 0 && manchester_values[i+1] == 0){ // LL (ignore)
        ind1++;
      }
    }

    for (int i = 1; i < manchester_values.size() - 1; i += 2){
      if (manchester_values[i] == 1 && manchester_values[i+1] == 1){ // HH (ignore)
        ind2++;
      }else if (manchester_values[i] == 0 && manchester_values[i+1] == 0){ // LL (ignore)
        ind2++;
      }
    }

    if (ind1 > ind2){
      initial = 1;
    }else{
      initial = 0;
    }
  }

	// in the case that we need state-saving we use this index
	// which will NOT be returned and end up overwriting the initial index
	// from block 1
	int start_index = initial;

	if (initial == 1 && !is_first_block){
		manchester_values.insert(manchester_values.begin(), cdr_state[0]);
		start_index = 0;
  }

	for (int i = start_index; i < manchester_values.size() - 1; i += 2){
		if (manchester_values[i] == 0 && manchester_values[i+1] == 1){	// LH = 0
			bits.push_back(0);
		}else if (manchester_values[i] == 1 && manchester_values[i+1] == 0){ // HL = 1
			bits.push_back(1);
		// change the two H's into one H
    }else if (manchester_values[i] == 1 && manchester_values[i+1] == 1){ // HH (ignore)
			//std::cerr << "Consecutive HI\n";
		//change the two L's into one L
    }else if (manchester_values[i] == 0 && manchester_values[i+1] == 0){ //LL (ignore)
			//std::cerr << "Consecutive LO\n";
    }
  }

	if (!is_first_block){
    bits.insert(bits.begin(), cdr_state[2]);
	}else{
		decoded_bits.push_back(bits[0]);
  }

	for (int i = 1; i < bits.size(); i++){
		decoded_bits.push_back(bits[i] ^ bits[i-1]);
  }

  return bits[bits.size() - 1];
}

std::tuple<char, int, std::vector<int>> frame_sync(std::vector<int> decoded_bits, std::vector<int> prev_decoded){
	int start_point = 0; //relative to bit stream
	if (prev_decoded.size() == 0){
		return {'X', start_point, decoded_bits};
	}else{
		//print("This is the obtained bit_stream")
		std::vector<int> bit_stream = prev_decoded;
    bit_stream.insert(bit_stream.end(), decoded_bits.begin(), decoded_bits.end());
  //print(bit_stream)
	//print(len(bit_stream))
  if (bit_stream.size() >=  26){
    std::vector<char> block = {'A','B','C','D'};
    for (int i = 0; i < bit_stream.size() - 25; i++){
    	std::vector<int> check = std::vector<int>(bit_stream.begin() + i, bit_stream.begin() + i + 26);
    	std::vector<int> message = matrix_mult(check, get_parity_check());
    	for (int k = 0; k < block.size(); k++){
    		if (message == get_syndrome(k)){
    			//std::cerr << "CORRECT SYNDROME OBTAINED\n";
    			//print(check)
    			//std::cerr << "obtained: " << block[k] << " " << i << "\n";
    			return {block[k], i, prev_decoded};
        }
      }
    }
  }
  return {'X', start_point, std::vector<int>(bit_stream.begin() + bit_stream.size() - 25, bit_stream.end())};
}
}

std::tuple<char, int, std::vector<int>> app_layer(char block_type, int start_point, std::vector<int> prev_decoded, std::vector<int> decoded_bits, std::string &d_service, int &d_index)
{
  // obtain bitstream to be parsed
  std::vector <int> bit_stream = prev_decoded;
  bit_stream.insert(bit_stream.end(), decoded_bits.begin(), decoded_bits.end());
  if (bit_stream.size() >= 26){
    //std::cerr << "RDS Block obtained\n";

    if (block_type == 'A'){
      //std::cerr << "---------------------BLOCK A---------------------\n";
      std::vector <int> pi_code = std::vector<int>(bit_stream.begin() + start_point, bit_stream.begin() + start_point + 16);
      std::string pi_hex = bin_to_hex(pi_code);
      //std::cerr << "The PI code is: " << pi_hex << "\n";

      block_type = 'B';
    }

    else if (block_type == 'B'){
      //std::cerr << "---------------------BLOCK B---------------------\n";
      std::vector <int> group_check = {0, 0, 0, 0, 0};
      std::vector <int> group_type = std::vector<int>(bit_stream.begin() + start_point, bit_stream.begin() + start_point + 5);
      std::vector <int> program_type;

      if (group_type == group_check)
      {
        program_type = std::vector<int>(bit_stream.begin() + start_point + 6, bit_stream.begin() + start_point + 11);
        std::cerr << "Group type is: ";
        for (int i = 0; i < group_type.size(); i++)
          std::cerr << group_type[i];

        std::cerr << "\nProgram type is ";
        for (int i = 0; i < program_type.size(); i++)
          std::cerr << program_type[i];
        std::cerr  << "\n";

        d_index = 2 * bit_stream[start_point + 14];
        d_index += bit_stream[start_point+15];
      }

      //else
        //std::cerr << "Group type is not 0A\n";

      block_type = 'C';
    }

    else if (block_type == 'C'){
      //std::cerr << "---------------------BLOCK C---------------------\n";
      // nothing happens

      block_type = 'D';
    }

    else if (block_type == 'D'){
      //std::cerr << "---------------------BLOCK D---------------------\n";
      std::vector <int> d1 = std::vector<int>(bit_stream.begin() + start_point, bit_stream.begin() + start_point + 8);
      std::vector <int> d2 = std::vector<int>(bit_stream.begin() + start_point + 8, bit_stream.begin() + start_point + 16);
      int d1_dec = bin_to_dec(d1);
      int d2_dec = bin_to_dec(d2);

      /*
      std::cerr << "D1 = ";
      for (int i = 0; i < d1.size(); i++)
        std::cerr << d1[i];
      std::cerr << "\n In decimal: " << d1_dec << "\n";
      std::cerr << "\n In ASCII: " << char(d1_dec) << "\n";

      std::cerr << "D2 = ";
      for (int i = 0; i < d2.size(); i++)
        std::cerr << d2[i];
      std::cerr << "\n In decimal: " << d2_dec << "\n";
      std::cerr << "\n In ASCII: " << char(d2_dec) << "\n";
      std::cerr << "d_index is " << d_index << "\n";
      */

      if (d_service.size() < 8)
      {
        if (d_index == 0 && d_service.size() == 0)
        {
          d_service.push_back(char(d1_dec));
          d_service.push_back(char(d2_dec));
        }

        else if (d_index == 1 && d_service.size() == 2)
        {
          d_service.push_back(char(d1_dec));
          d_service.push_back(char(d2_dec));
        }

        else if (d_index == 2 && d_service.size() == 4)
        {
          d_service.push_back(char(d1_dec));
          d_service.push_back(char(d2_dec));
        }

        else if (d_index == 3 && d_service.size() == 6)
        {
          d_service.push_back(char(d1_dec));
          d_service.push_back(char(d2_dec));
        }
      }

      else if (d_index == 0)
      {
        std::cerr << "Program Service: " << d_service << "\n";

        d_service.clear();
        d_service.push_back(char(d1_dec));
        d_service.push_back(char(d2_dec));
      }

      block_type = 'A';
    }

    else
      std::cerr << "Unknown block type!\n";

    return {block_type, 0, std::vector <int> (bit_stream.begin() + start_point + 26, bit_stream.end())};
  }

  else{
    start_point = 0;
    return {block_type, start_point, bit_stream};
  }
}

std::vector<int> matrix_mult(std::vector<int> x, std::vector<std::vector<int>> y){
	// This is not a generic matrix multiplication
	// This is build specificaly for bits and expected dimensions

	std::vector<int> result;
	int value = 0;
	// iterate through columns of y
	for (int i = 0; i < y[0].size(); i++){ // iterate through columns of x and rows of y

		for (int k = 0; k < x.size(); k++){
			//print(str(x[k]) + " and " + str(y[i][k]))
			//print(k)
			if (k == 0){
				value = (x[k] and y[k][i]);
			}else{
				value ^= (x[k] and y[k][i]);
      }
    }
		result.push_back(value);
  }
	return result;
}

std::vector<std::vector<int>> get_parity_check(){
  std::vector<std::vector<int>> parity_check = { \
  {1,0,0,0,0,0,0,0,0,0}, \
  {0,1,0,0,0,0,0,0,0,0}, \
  {0,0,1,0,0,0,0,0,0,0}, \
  {0,0,0,1,0,0,0,0,0,0}, \
  {0,0,0,0,1,0,0,0,0,0}, \
  {0,0,0,0,0,1,0,0,0,0}, \
  {0,0,0,0,0,0,1,0,0,0}, \
  {0,0,0,0,0,0,0,1,0,0}, \
  {0,0,0,0,0,0,0,0,1,0}, \
  {0,0,0,0,0,0,0,0,0,1}, \
  {1,0,1,1,0,1,1,1,0,0}, \
  {0,1,0,1,1,0,1,1,1,0}, \
  {0,0,1,0,1,1,0,1,1,1}, \
  {1,0,1,0,0,0,0,1,1,1}, \
  {1,1,1,0,0,1,1,1,1,1}, \
  {1,1,0,0,0,1,0,0,1,1}, \
  {1,1,0,1,0,1,0,1,0,1}, \
  {1,1,0,1,1,1,0,1,1,0}, \
  {0,1,1,0,1,1,1,0,1,1}, \
  {1,0,0,0,0,0,0,0,0,1}, \
  {1,1,1,1,0,1,1,1,0,0}, \
  {0,1,1,1,1,0,1,1,1,0}, \
  {0,0,1,1,1,1,0,1,1,1}, \
  {1,0,1,0,1,0,0,1,1,1}, \
  {1,1,1,0,0,0,1,1,1,1}, \
  {1,1,0,0,0,1,1,0,1,1}};
  return parity_check;
}

std::vector<int> get_syndrome(int row){
	// expected Syndrome array
	std::vector<std::vector<int>> exp_syndrome = { \
  {1,1,1,1,0,1,1,0,0,0}, \
	{1,1,1,1,0,1,0,1,0,0}, \
	{1,0,0,1,0,1,1,1,0,0}, \
	{1,0,0,1,0,1,1,0,0,0}};

	return exp_syndrome[row];
}

// built solely for block A, will not convert properly if bits not multiple of 4
std::string bin_to_hex(std::vector <int> bin)
{
    std::string hex = "";
    int bin_sum;

    for (int i = 0; i < bin.size(); i += 4)
    {
        bin_sum = 0;
        for (int k = 0; k < 4; k++)
            bin_sum += bin[i + k] * pow(2, 3 - k);

        if (bin_sum > 9)
        {
            if (bin_sum == 10)
                hex.push_back('A');
            else if (bin_sum == 11)
                hex.push_back('B');
            else if (bin_sum == 12)
                hex.push_back('C');
            else if (bin_sum == 13)
                hex.push_back('D');
            else if (bin_sum == 14)
                hex.push_back('E');
            else if (bin_sum == 15)
                hex.push_back('F');
        }

        else{
            char num = '0' + bin_sum;
            hex.push_back(num);
        }
    }
    return hex;
}

int bin_to_dec(std::vector <int> bin)
{
    int dec = 0;
    for (int i = bin.size() - 1; i >= 0; i--)
    {
        dec += bin[i]*pow(2, bin.size() - 1 - i);
    }
    return dec;
}
