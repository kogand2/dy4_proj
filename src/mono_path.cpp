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
void downsample(int D,std::vector<float> input, std::vector<float> &down)
{
    for(unsigned int i = 0; i < input.size(); i = i + D)
    {
        down.push_back(input[i]);
    }
}

std::vector<float> mono_path(std::vector<float> IQ_demod, std::vector<float> audio_coeff, std::vector<float> &audio_state, int audio_decim){
  // filter out audio data with convolution
  std::vector<float> audio_filt;
  state_block_conv(audio_filt, IQ_demod, audio_coeff, audio_state);

  // downsample filtered audio data
  std::vector<float> audio_block;
  downsample(audio_decim, audio_filt, audio_block);

  return audio_block;
}
