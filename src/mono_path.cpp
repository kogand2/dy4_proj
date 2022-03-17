
#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include "block_conv_fn.h"
#include "mono_path.h"

// read in raw block data
void readStdinBlockData(unsigned int num_samples, std::vector<float> &block_data){
  std::vector<char> raw_data(num_samples);
  std::cin.read(reinterpret_cast<char*>(&raw_data[0]), num_samples*sizeof(char));
  for(int k = 0; k < (int)num_samples; k++) {
    block_data[k] = float(((unsigned char)raw_data[k] - 128)/ 128.0);
  }
}

// dowsampling function
void downsample(int D, std::vector<float> input, std::vector<float> &down)
{
    // clear downsampled array
    down.clear();

    for(int i = 0; i < input.size(); i = i + D)
    {
        down.push_back(input[i]);
    }
}

std::vector<float> mono_path(int mode, std::vector<float> IQ_demod, std::vector<float> audio_coeff, std::vector<float> &audio_state, int audio_decim, int audio_exp){
  std::vector<float> audio_filt;
  std::vector<float> audio_block;

  if (mode == 2 || mode == 3)
  {
    //std::cerr << "test2\n";
    // resample audio data
    rs_block_conv(audio_filt, IQ_demod, audio_coeff, audio_state, audio_decim, audio_exp);


    // take downsampled filtered audio data
    audio_block.resize(audio_filt.size()/audio_decim);
    downsample(audio_decim, audio_filt, audio_block);
  }

  else{
    // filter out audio data with convolution
    ds_block_conv(audio_filt, IQ_demod, audio_coeff, audio_state, audio_decim);

    // take downsampled filtered audio data
    audio_block.resize(audio_filt.size()/audio_decim);
    downsample(audio_decim, audio_filt, audio_block);
  }

  return audio_block;
}
