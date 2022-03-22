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

std::vector<float> mono_path_i(int mode, std::vector<float> audio_filt, std::vector<float> IQ_demod, std::vector<float> audio_coeff, std::vector<float> &audio_state, int audio_decim, int audio_exp){
  std::vector<float> audio_block;

  if (mode == 2 || mode == 3)
  {
    // resample audio data
    //auto start_time = std::chrono::high_resolution_clock::now();
    audio_block.reserve(audio_filt.size()/audio_decim);
    rs_block_conv(audio_filt, IQ_demod, audio_coeff, audio_state, audio_decim, audio_exp, audio_block);

    // timing analysis
    /*auto stop_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> DFT_run_time = stop_time-start_time;
    std::cerr << "MONO RESAMPLE RUNTIME: " << DFT_run_time.count() << " ms" << "\n";
    */
    // take downsampled filtered audio data
    //downsample(audio_decim, audio_filt, audio_block);
  }

  else{
    // filter out audio data with convolution
    //auto start_time = std::chrono::high_resolution_clock::now();
    audio_block.reserve(audio_filt.size()/audio_decim);
    ds_block_conv(audio_filt, IQ_demod, audio_coeff, audio_state, audio_decim, audio_block);
    // timing analysis
    /*auto stop_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> DFT_run_time = stop_time-start_time;
    std::cerr << "MONO DOWNSAMPLE RUNTIME: " << DFT_run_time.count() << " ms" << "\n";
    */
  }

  return audio_block;
}
