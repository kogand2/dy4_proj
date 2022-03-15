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

// upsampling function
void upsample (int U, std::vector<float> input, std::vector<float> &up)
{
  // clear upsampled array
  up.clear();

  // pad with U - 1 zeros
  for(int i = 0; i < input.size(); i++)
  {
    for(int j = 0; j < U; j++)
      up.push_back(0);
  }
}

// convolution
void mono_convolution(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &state)
{
	// allocate memory for the output (filtered) data
	y.clear();
	y.resize(x.size(), 0.0);

  // lead-in
  for (int n = 0; n < h.size(); n++)
  {
    for (int k = 0; k < h.size(); k++){
        if (n-k >= 0)
          y[n] += h[k] * x[n-k];
        else
          y[n] += h[k] * state[state.size() + n - k];
    }
  }

  // dominant partition
  for (int n = h.size(); n < x.size(); n++)
  {
    for (int k = 0; k < h.size(); k++)
    {
        y[n] += h[k] * x[n-k];
    }
  }

  // lead out
  for (int n = x.size(); n < y.size(); n++)
  {
    for (int k = 0; k < h.size(); k++)
    {
      if (n-k < x.size())
        y[n] += h[k] * x[n-k];
      else
        y[n] += h[k] * state[state.size() + n - k];
    }
  }

  int index = x.size() - h.size() + 1;
	state = std::vector<float>(x.begin() + index, x.end());
}

std::vector<float> mono_path(int mode, std::vector<float> IQ_demod, std::vector<float> audio_coeff, std::vector<float> &audio_state, int audio_decim, int audio_exp){
  if (mode == 2 || mode == 3)
  {
    // upsample audio data
    std::vector<float> audio_us;
    audio_us.resize(IQ_demod.size()*audio_exp);
    upsample(audio_exp, IQ_demod, audio_us);

    // filter out audio data with convolution
    std::vector<float> audio_filt;
    state_block_conv(audio_filt, audio_us, audio_coeff, audio_state);

    // downsample filtered audio data
    std::vector<float> audio_block;
    audio_block.resize(audio_filt.size()/audio_decim);
    downsample(audio_decim, audio_filt, audio_block);

    return audio_block;
  }

  else{
  // filter out audio data with convolution
  std::vector<float> audio_filt;
  rf_block_conv(audio_filt, IQ_demod, audio_coeff, audio_state, audio_decim);

  // downsample filtered audio data
  std::vector<float> audio_block;
  audio_block.resize(audio_filt.size()/audio_decim);
  downsample(audio_decim, audio_filt, audio_block);

  return audio_block;
  }
}
