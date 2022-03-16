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
    // NOTE: THIS DESTROYS EVERYTHING IN VECTOR SO, THIS CAN BE FURTHER
    // OPTIMIZED BY RESIZING VECTOR ONCE AGAIN BEFORE PUSHING INTO IT
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
  up.resize(input.size()*U, 0.0);

  int k = 0; // used to parse through input

  // insert values in between each U - 1 zeros
  for(int i = 0; i < up.size(); i += U)
  {
    up[i] = input[k];
    k++;
  }
}

std::vector<float> mono_path(int mode, std::vector<float> IQ_demod, int audio_decim, int audio_exp, float audio_Fc, float audio_Fs, int audio_taps){
  if (mode == 2 || mode == 3)
  {
    // upsample audio data
    std::vector<float> audio_us;
    audio_us.resize(IQ_demod.size()*audio_exp, 0.0);
    upsample(audio_exp, IQ_demod, audio_us);

    // audio path LPF coefficients
  	std::vector<float> audio_coeff;
    int audio_Fs = (audio_exp*(audio_Fs/2))/audio_decim;
  	low_pass_coeff(audio_Fs, audio_Fc, audio_taps, audio_coeff);

    // state saving variable for audio data convolution
    std::vector<float> audio_state;
    audio_state.resize(audio_coeff.size(), 0.0);

    // filter out audio data with convolution
    std::vector<float> audio_filt;
    ds_block_conv(audio_filt, audio_us, audio_coeff, audio_state, audio_decim);

    // downsample filtered audio data
    std::vector<float> audio_block;
    audio_block.resize(audio_filt.size()/audio_decim);
    downsample(audio_decim, audio_filt, audio_block);

    return audio_block;
  }

  else{
    // no upsampling needed

    // audio path LPF coefficients
    std::vector<float> audio_coeff;
    low_pass_coeff(audio_Fs, audio_Fc, audio_taps, audio_coeff);

    // state saving variable for audio data convolution
  	std::vector<float> audio_state;
    audio_state.resize(audio_coeff.size(), 0.0);

    // filter out audio data with convolution
    std::vector<float> audio_filt;
    ds_block_conv(audio_filt, IQ_demod, audio_coeff, audio_state, audio_decim);

    // take downsampled filtered audio data
    std::vector<float> audio_block;
    audio_block.resize(audio_filt.size()/audio_decim);
    downsample(audio_decim, audio_filt, audio_block);

    return audio_block;
  }
}
