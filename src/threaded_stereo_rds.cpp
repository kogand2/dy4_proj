#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include "radio_fn.h"
#include "mono_path.h"

#include <stdio.h>
#include <math.h>

// testing time complexity
#include <chrono>
#include <iostream>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>

void readStdinBlockData(unsigned int num_samples, std::vector<float> &block_data){
  std::vector<char> raw_data(num_samples);
  std::cin.read(reinterpret_cast<char*>(&raw_data[0]), num_samples*sizeof(char));
  for(int k = 0; k < (int)num_samples; k++) {
    block_data[k] = float(((unsigned char)raw_data[k] - 128)/ 128.0);
  }
}

void back_end_rds_consumer(std::queue<std::vector<float>> &sync_queue, \
  std::mutex &mutex, \
  std::condition_variable &c_var,
  int &mode, \
  bool &end_of_stream)
{

  // Determine custom parameters depending on mode
  float rf_Fs;
  int block_size, rf_decim, sps, demod_decim, demod_exp;

  switch(mode) {
    case 0:
      rf_Fs = 2400000.0;
      rf_decim = 10;
      sps = 17;
      demod_decim = 1920;
      demod_exp = 323;
      block_size = 38400*2;
      break;
    case 2:
      rf_Fs = 2400000.0;
      rf_decim = 10;
      sps = 30;
      demod_decim = 64;
      demod_exp = 19;
      block_size = 24320*2;
      break;
  }
  // RF front end variables
  int audio_taps = 101;

  auto start_time = std::chrono::high_resolution_clock::now();
  auto stop_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> RDS_run_time;
  // RDS path variables
  int rds_sr = sps * 2375;

  // filter coefficients for RDS
  std::vector<float> carrier_coeff, channel_coeff, rds_demod_coeff, rds_rrc_coeff;

	band_pass_coeff(113500, 114500, rf_Fs/rf_decim, audio_taps, carrier_coeff);
	band_pass_coeff(54000, 60000, rf_Fs/rf_decim, audio_taps, channel_coeff);
	low_pass_coeff((rf_Fs/rf_decim)*demod_exp, 3000, audio_taps*demod_exp, rds_demod_coeff);
  rrc_coeff(rds_sr, audio_taps, rds_rrc_coeff);

  // RDS filtered data
  std::vector<float> carrier_filt, channel_filt, recoveredRDS, rds_demod_filt, rrc_demod_filt, demod_down, rds_demod_block;
  carrier_filt.resize((block_size / 2) - 1);
  channel_filt.resize((block_size / 2) - 1);
  rds_demod_filt.resize(demod_exp*((block_size / 2) - 1));
  rds_demod_block.resize(rds_demod_filt.size()/demod_decim);

  rrc_demod_filt.resize(rds_demod_block.size());

  // state saving variables for RDS path
	std::vector<float> carrier_state, channel_state, rds_demod_state, rds_rrc_state, pll_state;
	carrier_state.resize(carrier_coeff.size() - 1);	// for RDS Carrier Recovery
	channel_state.resize(channel_coeff.size() - 1);	// for RDS Channel Extraction
  rds_demod_state.resize(rds_demod_coeff.size() - 1); // for RDS Rational Resampler
  rds_rrc_state.resize(rds_rrc_coeff.size() - 1); // for RDS RRC convolution

  pll_state.resize(6);
  pll_state[0] = 0.0;
  pll_state[1] = 0.0;
  pll_state[2] = 1.0;
  pll_state[3] = 0.0;
  pll_state[4] = 1.0;
  pll_state[5] = 0.0;

  // RDS demodulation
  std::vector<float> mixed_audio, delay_block, channel_delay;
  mixed_audio.resize((block_size / 2) - 1);
  delay_block.resize((audio_taps - 1)/2);

  float prev_phase = 0.0;

  int cdr_init = -1;
  int decode_init = -1;
  std::vector<int> cdr_state = {0, 0, 0};
  std::vector<int> sample_idxs, man_encoding, decoded_bits, prev_decoded;
  int last_bit;

  char block_type = 'X';
  int start_point = 0;
  std::string d_service;
  int d_index = 0;
  bool is_first_block = true;

  while (!end_of_stream || !sync_queue.empty()){
    // STEP 2: Stereo path
    std::unique_lock<std::mutex> lock(mutex);

    while (sync_queue.empty()){
      c_var.wait(lock);
    }

    std::vector<float> IQ_demod = sync_queue.front();

    sync_queue.pop();

    if (mode == 0 || mode == 2){
      c_var.notify_all();
    }else{
      c_var.notify_one();
    }

    lock.unlock();

    // STEP 2: RDS path
    start_time = std::chrono::high_resolution_clock::now();

    //RDS Channel Extraction: all-pass
		state_block_conv(channel_filt, IQ_demod, channel_coeff, channel_state);
    all_pass_coeff(channel_delay, channel_filt, delay_block);

		//RDS Carrier Recovery: non-linearity -> Bandpass -> PLL -> Numerically Controlled
    std::vector<float> carrier_input = sq_non_linearity(channel_filt);

		state_block_conv(carrier_filt, carrier_input, carrier_coeff, carrier_state);

		fmPll(carrier_filt, recoveredRDS, pll_state, 114e3, rf_Fs/rf_decim, 0.5, 1.1775, 0.003);

    //RDS Demodulation: mixer -> digital filtering (lowpass -> resampler) -> rrc -> CDR
		mixer(recoveredRDS, channel_delay, mixed_audio);

    rds_demod_block.clear();
    rs_block_conv(rds_demod_filt, mixed_audio, rds_demod_coeff, rds_demod_state, demod_decim, demod_exp, rds_demod_block);
    state_block_conv(rrc_demod_filt, rds_demod_block, rds_rrc_coeff, rds_rrc_state);

    // RDS Demodulation and Data Processing: CDR, decoding, frame snychronization, app layer
    if (!is_first_block){
      CDR(rrc_demod_filt, sps, cdr_init, sample_idxs, man_encoding);
      last_bit = diff_decoding(man_encoding, cdr_state, decoded_bits, decode_init);
      //printRealVector(decoded_bits);
      cdr_state = {man_encoding[man_encoding.size() - 1], sample_idxs[sample_idxs.size() - 1], last_bit};

      if (block_type == 'X'){
        std::tie(block_type, start_point, prev_decoded) = frame_sync(decoded_bits, prev_decoded);
        // frame synchronized: can begin information extraction

        if (block_type != 'X')
          std::tie(block_type, start_point, prev_decoded) = app_layer(block_type, start_point, prev_decoded, decoded_bits, d_service, d_index);
      }

      else
        std::tie(block_type, start_point, prev_decoded) = app_layer(block_type, start_point, prev_decoded, decoded_bits, d_service, d_index);
    }
    is_first_block = false;

    // timing analysis
    stop_time = std::chrono::high_resolution_clock::now();
    RDS_run_time += stop_time-start_time;
  }
  std::cerr << "RDS PATH RUNTIME: " << RDS_run_time.count() << " ms" << "\n";
}

void back_end_stereo_consumer(std::queue<std::vector<float>> &sync_queue, \
  std::mutex &mutex, \
  std::condition_variable &c_var,
  int &mode, \
  bool &end_of_stream)
{

  // Determine custom parameters depending on mode
  float rf_Fs;
  int block_size, rf_decim, audio_decim, audio_exp;

  switch(mode) {
    case 0:
      rf_Fs = 2400000.0;
      rf_decim = 10;
      audio_exp = 1;
      audio_decim = 5;
      block_size = 38400*2;
      break;
    case 1:
      rf_Fs = 1440000.0;
      rf_decim = 5;
      audio_exp = 1;
      audio_decim = 6;
      block_size = 102400; //needs to be changed maybe
      break;
    case 2:
      rf_Fs = 2400000.0;
      rf_decim = 10;
      audio_exp = 147;
      audio_decim = 800;
      block_size = 24320*2;
      break;
    case 3:
      rf_Fs = 2304000.0;
      rf_decim = 9;
      audio_exp = 441;
      audio_decim = 2560;
      block_size = 38400*2; //needs to be changed maybe
      break;
  }

  // audio path variables
  float audio_Fc = 16000;
  int audio_taps = 101;

  auto start_time = std::chrono::high_resolution_clock::now();
  auto stop_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> STEREO_run_time;

  std::vector<float> carrier_coeff, channel_coeff, audio_coeff;
  band_pass_coeff(18500, 19500, rf_Fs/rf_decim, audio_taps, carrier_coeff);
  band_pass_coeff(22000, 54000, rf_Fs/rf_decim, audio_taps, channel_coeff);
  low_pass_coeff((rf_Fs/rf_decim)*audio_exp, audio_Fc, audio_taps*audio_exp, audio_coeff);

  // stereo carrier recovery, channel extraction, and mixed audio
  std::vector<float> carrier_filt, channel_filt, recoveredStereo;
  carrier_filt.resize((block_size / 2) - 1);
  channel_filt.resize((block_size / 2) - 1);

  // state saving variables for stereo path
	std::vector<float> stereo_state, carrier_state, channel_state, pll_state;
	stereo_state.resize(carrier_coeff.size() - 1);	// for Stereo Processing
	carrier_state.resize(carrier_coeff.size() - 1);	// for Stereo Carrier Recovery
	channel_state.resize(channel_coeff.size() - 1);	// for Stereo Channel Extraction

  pll_state.resize(6);
  pll_state[0] = 0.0;
  pll_state[1] = 0.0;
  pll_state[2] = 1.0;
  pll_state[3] = 0.0;
  pll_state[4] = 1.0;
  pll_state[5] = 0.0;

  // stereo processing
  std::vector<float> mixed_audio, stereo_filt, stereo_block;
  mixed_audio.resize((block_size / 2) - 1);
  stereo_filt.resize((mixed_audio.size()-1)*audio_exp);
  stereo_block.resize(stereo_filt.size()/audio_decim);

  // left and right audio data
  float left_block, right_block;
	std::vector<float> delay_block, mono_input;
  delay_block.resize((audio_taps - 1)/2);

  std::vector<std::vector<float>> complete_data;

  // filtered audio data
  std::vector<float> audio_filt;
  audio_filt.resize(audio_exp*(((block_size / 2) / rf_decim) - 1));

  // downsampled filtered audio data
  std::vector<float> audio_block;
  audio_block.resize(audio_filt.size()/audio_decim);

  // state saving variable for audio data convolution
  std::vector<float> audio_state;
  audio_state.resize(audio_coeff.size() - 1);
  int count = 0;
  while (!end_of_stream || !sync_queue.empty()){
    // STEP 2: Stereo path
    std::unique_lock<std::mutex> lock(mutex);

    while (sync_queue.empty()){
      c_var.wait(lock);
    }

    std::vector<float> IQ_demod = sync_queue.front();

    sync_queue.pop();

    if (mode == 0 || mode == 2){
      c_var.notify_all();
    }else{
      c_var.notify_one();
    }

    lock.unlock();
    count += 1;

    std::cerr << "count: " << count << "\n";
    start_time = std::chrono::high_resolution_clock::now();
    //Stereo Carrier Recovery: Bandpass -> PLL -> Numerically Controlled
		state_block_conv(carrier_filt, IQ_demod, carrier_coeff, carrier_state);
		fmPll(carrier_filt, recoveredStereo, pll_state, 19e3, rf_Fs/rf_decim, 2.0, 0.0, 0.01);

    //Stereo Channel Extraction: Bandpass
		state_block_conv(channel_filt, IQ_demod, channel_coeff, channel_state);

    //Stereo Processing: Mixer -> Digital filtering (Lowpass -> down sample) -> Stereo Combiner
		mixer(recoveredStereo, channel_filt, mixed_audio);

    stereo_block.clear();
    std::cerr << "test1\n";
    if (mode == 0 || mode == 1)
		  ds_block_conv(stereo_filt, mixed_audio, audio_coeff, stereo_state, audio_decim, stereo_block);
    else
      rs_block_conv(stereo_filt, mixed_audio, audio_coeff, stereo_state, audio_decim, audio_exp, stereo_block);

    all_pass_coeff(mono_input, IQ_demod, delay_block);
    std::cerr << "test2\n";
    audio_block.clear();
    if (mode == 0 || mode == 1)
      ds_block_conv(audio_filt, mono_input, audio_coeff, audio_state, audio_decim, audio_block);
    else
      rs_block_conv(audio_filt, mono_input, audio_coeff, audio_state, audio_decim, audio_exp, audio_block);

    // timing analysis
    stop_time = std::chrono::high_resolution_clock::now();
    STEREO_run_time += stop_time-start_time;

    // STEP 3: prepare audio data for output
    std::vector<short int> audio_data;

    for (unsigned int k = 0; k < audio_block.size(); k++) {
      left_block = audio_block[k] + stereo_block[k];
			right_block = audio_block[k] - stereo_block[k];

      if(std::isnan(left_block)){
        audio_data.push_back(0);
        audio_data.push_back(0);
      }
      else{
				audio_data.push_back(static_cast<short int>(left_block * 16384));
				audio_data.push_back(static_cast<short int>(right_block * 16384));
      }
    }

    fwrite(&audio_data[0], sizeof(short int), audio_data.size(), stdout);
  }
  std::cerr << "test3\n";
  std::cerr << "STEREO PATH RUNTIME: " << STEREO_run_time.count() << " ms" << "\n";
}

void rf_front_end_producer(std::queue<std::vector<float>> &sync_queue_stereo, \
  std::queue<std::vector<float>> &sync_queue_rds, \
  std::mutex &mutex, \
  std::condition_variable &c_var,
  int &mode, \
  bool &end_of_stream)
{
  // TIMING VARIABLES
  auto start_time = std::chrono::high_resolution_clock::now();
  auto stop_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> RF_run_time;

  // Determine custom parameters depending on mode
  float rf_Fs, audio_Fs;
  int block_size, rf_decim, audio_decim, audio_exp;

  switch(mode) {
    case 0:
      rf_Fs = 2400000.0;
      rf_decim = 10;
      block_size = 38400*2;
      break;
    case 1:
      rf_Fs = 1440000.0;
      rf_decim = 5;
      block_size = 102400; //needs to be changed maybe?
      break;
    case 2:
      rf_Fs = 2400000.0;
      rf_decim = 10;
      block_size = 24320*2;
      break;
    case 3:
      rf_Fs = 2304000.0;
      rf_decim = 9;
      block_size = 38400*2; //needs to be changed maybe?
      break;
  }

  // RF front end variables
	float rf_Fc = 100000.0;
	int rf_taps = 101;

  // filter coefficients for RF
  std::vector<float> rf_coeff;
	low_pass_coeff(rf_Fs, rf_Fc, rf_taps, rf_coeff);
  std::cerr << "test\n";
  // regular IQ data
  std::vector<float> i_data, q_data;
  i_data.resize(block_size / 2);
  q_data.resize(block_size / 2);

  // filtered IQ data
  std::vector<float> i_filt, q_filt;
  i_filt.resize(i_data.size());
  q_filt.resize(q_data.size());

  // downsampled filtered IQ data
  std::vector<float>i_ds, q_ds;
  i_ds.resize(i_filt.size()/rf_decim);
  q_ds.resize(q_filt.size()/rf_decim);

  // state saving variables for I and Q samples convolution
  std::vector<float> state_i_lpf_100k, state_q_lpf_100k;
  state_i_lpf_100k.resize(rf_coeff.size() - 1, 0.0);
  state_q_lpf_100k.resize(rf_coeff.size() - 1, 0.0);

  // demodulation variables
	std::vector<float> IQ_demod;
  float prev_phase = 0.0;
  int block_id = 0;

  // decipher each block
	while(1) {
    //std::cerr << "Produce\n";
		std::vector<float> block_data(block_size);
    readStdinBlockData(block_size, block_data);
    if ((std::cin.rdstate()) != 0) {
      end_of_stream = true;
      std::cerr << "End of input stream reached" << std::endl;
      std::cerr << "RF DOWNSAMPLE RUNTIME: " << RF_run_time.count() << " ms" << "\n";
      break;
    }

    // STEP 1: IQ samples demodulation
    // take IQ data in
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

    // perform demodulation on IQ data
    IQ_demod = fmDemodArctan(i_ds, q_ds, prev_phase);

    // timing analysis
    stop_time = std::chrono::high_resolution_clock::now();
    RF_run_time += stop_time-start_time;

    std::unique_lock<std::mutex> lock(mutex);

    if (mode == 0 || mode == 2){
      while (sync_queue_stereo.size() >= 5 || sync_queue_rds.size() >= 5){
        c_var.wait(lock);
      }

      sync_queue_stereo.push(IQ_demod);
      sync_queue_rds.push(IQ_demod);

      c_var.notify_all();
      lock.unlock();
    }else{
      while (sync_queue_stereo.size() >= 5){
        c_var.wait(lock);
      }

      sync_queue_stereo.push(IQ_demod);

      c_var.notify_one();
      lock.unlock();
    }
    block_id++;
  }
}

void stereo_path_main(int mode){
	// TIMING VARIABLES
  auto start_time = std::chrono::high_resolution_clock::now();
  bool end_of_stream = false;
  std::queue<std::vector<float>> sync_queue_stereo, sync_queue_rds;
  std::mutex mutex;
  std::condition_variable c_var;

  if (mode == 0 || mode == 2){
    std::thread rf_front_end = std::thread(rf_front_end_producer, std::ref(sync_queue_stereo), \
      std::ref(sync_queue_rds), std::ref(mutex), std::ref(c_var), std::ref(mode), std::ref(end_of_stream));

    std::thread stereo_back_end = std::thread(back_end_stereo_consumer, std::ref(sync_queue_stereo), \
      std::ref(mutex), std::ref(c_var), std::ref(mode), std::ref(end_of_stream));

    std::thread rds_back_end = std::thread(back_end_rds_consumer, std::ref(sync_queue_rds), \
      std::ref(mutex), std::ref(c_var), std::ref(mode), std::ref(end_of_stream));

    rf_front_end.join();
    stereo_back_end.join();
    rds_back_end.join();

  }else{
    std::cerr << "RDS IS NOT AVAILABLE IN THIS MODE.\n";

    std::thread rf_front_end = std::thread(rf_front_end_producer, std::ref(sync_queue_stereo), \
      std::ref(sync_queue_rds), std::ref(mutex), std::ref(c_var), std::ref(mode), std::ref(end_of_stream));

    std::thread stereo_back_end = std::thread(back_end_stereo_consumer, std::ref(sync_queue_stereo), \
      std::ref(mutex), std::ref(c_var), std::ref(mode), std::ref(end_of_stream));

    rf_front_end.join();
    stereo_back_end.join();
  }

  auto stop_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> total_time = stop_time - start_time;
  std::cerr << "TOTAL TIME: " << total_time.count() << " ms" << "\n";
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
