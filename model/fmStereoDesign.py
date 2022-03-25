#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy import signal
import numpy as np
import math

# use fmDemodArctan and fmPlotPSD
from fmSupportLib import *
#from fmMonoBlock import *
# for take-home add your functions


rf_Fs = 2.4e6
rf_Fc = 100e3
rf_taps = 151
rf_decim = 10

audio_Fs = 48e3
audio_decim = 5
audio_taps = 151
audio_Fc = 16e3
# add other settings for audio, like filter taps, ...

# flag that keeps track if your code is running for
# in-lab (il_vs_th = 0) vs takehome (il_vs_th = 1)

if __name__ == "__main__":

	# read the raw IQ data from the recorded file
	# IQ data is assumed to be in 8-bits unsigned (and interleaved)
	in_fname = "../data/samples0.raw"
	raw_data = np.fromfile(in_fname, dtype='uint8')
	print("Read raw RF data from \"" + in_fname + "\" in unsigned 8-bit format")
	# IQ data is normalized between -1 and +1 in 32-bit float format
	iq_data = (np.float32(raw_data) - 128.0)/128.0
	print("Reformatted raw RF data to 32-bit float format (" + str(iq_data.size * iq_data.itemsize) + " bytes)")

	# coefficients for the front-end low-pass filter
	rf_coeff = signal.firwin(rf_taps, rf_Fc/(rf_Fs/2), window=('hann'))

	carrier_coeff = myBandPass(18.5e3, 19.5e3, rf_Fs/rf_decim, audio_taps) #For Stereo Carrier Recovery

	#th_coeff = signal.firwin(audio_taps, [18.5e3, 19.5e3], pass_zero=False, fs = rf_Fs/rf_decim) #For testing

	#carrier_coeff1 = myLowPass(18.5e3, rf_Fs/rf_decim, audio_taps) #For mono path combination #For testing
	#carrier_coeff2 = myLowPass(19.5e3, rf_Fs/rf_decim, audio_taps) #For mono path combination
	#Low_Carrier_Coeff = np.subtract(carrier_coeff2,carrier_coeff1)

	#x1 = range(151)
	#fig2, axs = plt.subplots(3)
	#fig2.suptitle('Coeff Checking')
	#axs[0].plot(x1, carrier_coeff)
	#axs[1].plot(x1, th_coeff)
	#axs[2].plot(x1, Low_Carrier_Coeff)
	#plt.show()

	channel_coeff = myBandPass(22e3, 54e3, rf_Fs/rf_decim, audio_taps) #For Stereo Channel Extraction

	audio_coeff = myLowPass(audio_Fc, rf_Fs/rf_decim, audio_taps) #For mono path combination
	mono_coeff = myLowPass(audio_Fc, rf_Fs/rf_decim, audio_taps) #For mono path combination

	# set up the subfigures for plotting
	subfig_height = np.array([0.8, 2, 2, 1.6]) # relative heights of the subfigures
	plt.rc('figure', figsize=(10, 10))	# the size of the entire figure
	fig, (ax0, ax1, ax3, ax2) = plt.subplots(nrows=4, gridspec_kw={'height_ratios': subfig_height})
	fig.subplots_adjust(hspace = .6)

	# select a block_size that is a multiple of KB
	# and a multiple of decimation factors
	block_size = 1024 * rf_decim * audio_decim * 2
	block_count = 0

	# states needed for continuity in block processing
	state_i_lpf_100k = np.zeros(rf_taps-1)
	state_q_lpf_100k = np.zeros(rf_taps-1)
	state_phase = 0

	dummy_state = np.array([0,0])						#For RF Demondulation
	filt_block = np.zeros(shape=len(audio_coeff))		#For Stereo Processing
	carrier_block = np.zeros(shape=len(carrier_coeff))	#For Stereo Carrier Recovery
	channel_block = np.zeros(shape=len(channel_coeff))	#For Stereo Channel Extraction
	mono_block = np.zeros(shape=len(mono_coeff))		#For Monopath

	left_block = np.zeros(shape=1024)		#For Monopath
	right_block = np.zeros(shape=1024)		#For Monopath
	delay_block = np.zeros(shape=(audio_taps-1)//2)		#For Monopath


	pll_block = np.array([0.0, 0.0, 1.0, 0.0, 1.0, 0.0])
	# add state as needed for the mono channel filter

	prevStereo = np.zeros(shape=5121)
	prevCarrier = np.zeros(shape=5120)

	prevMono = np.zeros(shape=5120)
	prevFilter = np.zeros(shape=5120)

	# audio buffer that stores all the audio blocks
	audio_data = np.array([]) # used to concatenate filtered blocks (audio data)
	stereo_data = np.array([])# used to concatenate filtered blocks (stereo data)

	left_data = np.array([]) # used to concatenate filtered left blocks (audio data)
	right_data = np.array([])# used to concatenate filtered right blocks (stereo data)

	while (block_count+1)*block_size < len(iq_data):
	#RF FRONT END==================================================
		# if you wish to have shorter runtimes while troubleshooting
		# you can control the above loop exit condition as you see fit
		print('Processing block ' + str(block_count))

		# filter to extract the FM channel (I samples are even, Q samples are odd)
		i_filt, state_i_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
				iq_data[(block_count)*block_size:(block_count+1)*block_size:2],
				zi=state_i_lpf_100k)
		q_filt, state_q_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
				iq_data[(block_count)*block_size+1:(block_count+1)*block_size:2],
				zi=state_q_lpf_100k)

		# downsample the I/Q data from the FM channel
		i_ds = i_filt[::rf_decim]
		q_ds = q_filt[::rf_decim]

		# you MUST have also "custom" state-saving for your own FM demodulator
		dummy_fm, state_phase = fmDemodArctan(i_ds, q_ds, prev_phase = state_phase)
		#dummy_fm, dummy_state = myDemod(i_ds, q_ds, dummy_state)

		#print(len(i_ds))
    #RF front end stops here==================================================

		# extract the stereo audio data

        #Stereo Carrier Recovery: Bandpass -> PLL -> Numerically Controlled Oscillator
		carrier_filt, carrier_block = block_convolution(carrier_coeff, dummy_fm, carrier_block)


		recoveredStereo, pll_block = fmPll(carrier_filt, 19e3, rf_Fs/rf_decim, 2.0, 0.0, 0.01, pll_block)


        #Stereo Channel Extraction: Bandpass
		channel_filt,channel_block = block_convolution(channel_coeff, dummy_fm, channel_block)


        #Stereo Processing: Mixer -> Digital filtering (Lowpass -> down sample) -> Stereo Combiner
		mixedAudio = mixer(recoveredStereo, channel_filt)

		stereo_block,filt_block = ds_block_convolution(audio_coeff, mixedAudio, filt_block, audio_decim)

		# extract the mono audio data
		mono_input, delay_block = myAllPass(dummy_fm, delay_block)

		audio_block,mono_block = ds_block_convolution(audio_coeff, mono_input, mono_block, audio_decim)


		for i in range(len(audio_block)):
			left_block[i] = (audio_block[i] + stereo_block[i])
			#print("Left Block [" + str(i) + "] = " + str(left_block[i]*4095))
			right_block[i] = (audio_block[i] - stereo_block[i])
			#print("Right Block [" + str(i) + "] = " + str(right_block[i]*4095))



		left_data = np.concatenate((left_data, left_block))
		right_data = np.concatenate((right_data, right_block))

		#Generate Plots of PLL
		#if block_count >= 3 and block_count < 6:
		#	x1 = range(50)
		#	x2 = range(50)
		#	fig2, axs = plt.subplots(2)
		#	fig2.suptitle('PLL Checking')
		#	axs[0].plot(x2, carrier_filt[:50], c='blue')
		#	axs[0].plot(x2, prevCarrier[5070:], c='orange')
		#	axs[0].set_title('Carrier Input', fontstyle='italic',fontsize='medium')
		#	axs[1].plot(x1, recoveredStereo[:50], c='blue')
		############################	axs[1].plot(x1, prevStereo[5071:], c='orange')
		#	axs[1].set_title('Carrier Output', fontstyle='italic',fontsize='medium')
		#	plt.show()
		#prevStereo = recoveredStereo
		#prevCarrier = carrier_filt

		#Generate Plots of Monopath
		if block_count >= 3 and block_count < 6:
			#print("hello")
			print(len(prevFilter[974:]))
			x1 = range(50)
			x2 = range(50)
			fig2, axs = plt.subplots(3)
			fig2.suptitle('State saving checking')
			axs[0].plot(range(50,100), stereo_block[:50], c='blue')
			axs[0].plot(x2, prevFilter[974:], c='orange')
			axs[0].set_title('Stereo path', fontstyle='italic',fontsize='medium')

			axs[1].plot(range(50,100), mono_block[:50], c='blue')
			axs[1].plot(x1, prevMono[(len(prevMono)-50):], c='orange')
			axs[1].set_title('Mono path', fontstyle='italic',fontsize='medium')

			axs[2].plot(range(49,99), recoveredStereo[:50], c='blue')
			axs[2].plot(x1, prevStereo[5071:], c='orange')
			axs[2].set_title('PLL', fontstyle='italic',fontsize='medium')

			plt.show()

		prevMono = mono_block
		prevFilter = stereo_block
		prevStereo = recoveredStereo

		if block_count >= 10 and block_count < 12:

			# plot PSD of selected block after FM demodulation
			ax0.clear()
			fmPlotPSD(ax0, dummy_fm, (rf_Fs/rf_decim)/1e3, subfig_height[0],
					'Demodulated FM (block ' + str(block_count) + ')')
			# output binary file name (where samples are written from Python)
			fm_demod_fname = "../data/fm_demod_" + str(block_count) + ".bin"
			# create binary file where each sample is a 32-bit float
			dummy_fm.astype('float32').tofile(fm_demod_fname)

			# plot PSD of selected block after extracting mono audio
			# ... change as needed
			fmPlotPSD(ax1, audio_block, ((rf_Fs/rf_decim)/1e3), subfig_height[1], 'mixed audio')

			fmPlotPSD(ax2, right_data, ((rf_Fs/rf_decim)/1e3), subfig_height[2], 'Right Channel')
			# plot PSD of selected block after downsampling mono audio
			# ... change as needed
			fmPlotPSD(ax3, left_data, audio_Fs/1e3, subfig_height[3], 'Left Channel')
			# save figure to file
			fig.savefig("../data/fmMonoBlockDelay" + str(block_count) + ".png")

		block_count += 1

	print('Finished processing all the blocks from the recorded I/Q samples')

	complete_data = np.vstack((left_data, right_data)).T

	# write audio data to file
	out_fname = "../data/fmStereoOldDemod.wav"
	wavfile.write(out_fname, int(audio_Fs), np.int16((complete_data/2)*4095))
	print("Written audio samples to \"" + out_fname + "\" in signed 16-bit format")

	# uncomment assuming you wish to show some plots
	# plt.show()
