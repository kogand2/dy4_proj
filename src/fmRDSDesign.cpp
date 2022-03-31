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

std::vector<float> rds_path_main(int mode){
	// TIMING VARIABLES
  auto start_overall_time = std::chrono::high_resolution_clock::now();
  auto start_time = std::chrono::high_resolution_clock::now();
  auto stop_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> RDS_run_time;
  std::chrono::duration<double, std::milli> RF_run_time;

  // Determine custom parameters depending on mode
  float rf_Fs;
  int rf_decim, sps, demod_decim, demod_exp;

  switch(mode) {
    case 0:
      rf_Fs = 2400000.0;
      rf_decim = 10;
      sps = 17;
      demod_decim = 1920;
      demod_exp = 323;
      break;
    case 1:
      rf_Fs = 1440000.0;
      rf_decim = 5;
      // HANDLE WITH PRINT STATEMENT
      break;
    case 2:
      rf_Fs = 2400000.0;
      rf_decim = 10;
      sps = 30;
      demod_decim = 1; // NOT CALCULATED
      demod_exp = 1; // NOT CALCULATED
      break;
    case 3:
      rf_Fs = 2304000.0;
      rf_decim = 9;
      // HANDLE WITH PRINT STATEMENT
      break;
  }

  // RF front end variables
	float rf_Fc = 100000.0;
	int rf_taps = 151;
  int audio_taps = 151;
  // RDS path variables
  int rds_sr = sps * 2375;

  // filter coefficients for RDS
  std::vector<float> rf_coeff, carrier_coeff, channel_coeff, rds_demod_coeff, rds_rrc_coeff;
	low_pass_coeff(rf_Fs, rf_Fc, rf_taps, rf_coeff);
	band_pass_coeff(113500, 114500, rf_Fs/rf_decim, audio_taps, carrier_coeff);
	band_pass_coeff(54000, 60000, rf_Fs/rf_decim, audio_taps, channel_coeff);
	low_pass_coeff((rf_Fs/rf_decim)*demod_exp, 3000, audio_taps*demod_exp, rds_demod_coeff);
  rrc_coeff(rds_sr, rf_taps, rds_rrc_coeff);
	int block_size = 38400*2;

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
	std::vector<float> demod_state, IQ_demod;
	demod_state.resize(2, 0.0);

  // RDS filtered data
  std::vector<float> carrier_filt, channel_filt, recoveredRDS, rds_demod_filt, rrc_demod_filt, demod_down, rds_demod_block;
  carrier_filt.resize(i_data.size() - 1);
  channel_filt.resize(i_data.size() - 1);
  rds_demod_filt.resize(demod_exp*(i_data.size() - 1));
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
  mixed_audio.resize(i_data.size() - 1);
  delay_block.resize((rf_taps - 1)/2);

  float prev_phase = 0.0;

  int cdr_init = -1;
  int decode_init = -1;
  std::vector<int> cdr_state = {0, 0, 0};
  std::vector<int> sample_idxs, man_encoding, decoded_bits;
  int last_bit;
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
      std::cerr << "TOTAL RDSPATH RUNTIME: " << RDS_run_time.count() << " ms" << "\n";
      exit(1);
    }

    //std::cerr << "Processing Block " << block_id << std::endl;

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

    if (block_id >= 1){
      CDR(rrc_demod_filt, sps, cdr_init, sample_idxs, man_encoding);
      last_bit = diff_decoding(man_encoding, cdr_state, decoded_bits, decode_init);
      printRealVector(decoded_bits);
      cdr_state = {man_encoding[man_encoding.size() - 1], sample_idxs[sample_idxs.size() - 1], last_bit};
    }
    //...

    // timing analysis
    stop_time = std::chrono::high_resolution_clock::now();
    RDS_run_time += stop_time-start_time;
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
  rds_path_main(mode);

	return 0;
}
