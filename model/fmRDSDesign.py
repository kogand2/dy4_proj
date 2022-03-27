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
audio_Fc = 3e3

# RDS VARIABLES=================================================================
# GCD = 125, U = 2375*17/125, D = 240000*U/(2375*17)
sps = 17
rds_sr = sps * 2375
demod_decim = 1920
demod_exp = 323

# add other settings for audio, like filter taps, ...

if __name__ == "__main__":

	# read the raw IQ data from the recorded file
	# IQ data is assumed to be in 8-bits unsigned (and interleaved)
	in_fname = "../data/samples0_2400.raw"
	raw_data = np.fromfile(in_fname, dtype='uint8')
	print("Read raw RF data from \"" + in_fname + "\" in unsigned 8-bit format")
	# IQ data is normalized between -1 and +1 in 32-bit float format
	iq_data = (np.float32(raw_data) - 128.0)/128.0
	print("Reformatted raw RF data to 32-bit float format (" + str(iq_data.size * iq_data.itemsize) + " bytes)")

	# set up the subfigures for plotting
	subfig_height = np.array([0.8, 2, 2, 1.6]) # relative heights of the subfigures
	plt.rc('figure', figsize=(10, 10))	# the size of the entire figure
	fig, (ax0, ax1, ax3, ax2) = plt.subplots(nrows=4, gridspec_kw={'height_ratios': subfig_height})
	fig.subplots_adjust(hspace = .6)

	# FILTER COEFFICIENTS
	rf_coeff = signal.firwin(rf_taps, rf_Fc/(rf_Fs/2), window=('hann'))
	carrier_coeff = myBandPass(113.5e3, 114.5e3, rf_Fs/rf_decim, audio_taps) #For RDS Carrier Recovery
	channel_coeff = myBandPass(54e3, 60e3, rf_Fs/rf_decim, audio_taps) #For RDS Channel Extraction
	rds_demod_coeff = myLowPass(3e3, rf_Fs/rf_decim, audio_taps)	# For RDS demodulation
	rrc_coeff = rrc(rds_sr, audio_taps) # for RDS demodulation

	# select a block_size that is a multiple of KB
	# and a multiple of decimation factors
	block_size = 1024 * rf_decim * audio_decim * 2
	block_count = 0

	# STATE SAVING IQ
	state_i_lpf_100k = np.zeros(rf_taps-1)
	state_q_lpf_100k = np.zeros(rf_taps-1)
	state_phase = 0

	dummy_state = np.array([0,0])						#For RF Demondulation

	# FILTERED DATA
	channel_filt = np.zeros(shape=block_size//2 - 1)
	carrier_filt = np.zeros(shape=block_size//2 - 1)
	channel_Delay = np.zeros(shape=block_size//2 - 1)
	demod_filt = np.zeros(shape=block_size//2 - 1)

	# STATE SAVING
	carrier_block = np.zeros(shape=len(carrier_coeff))	#For Carrier Recovery
	channel_block = np.zeros(shape=len(channel_coeff))	#For Channel Extraction
	delay_block = np.zeros(shape=((audio_taps-1)//2))		#For All pass filter
	rds_demod_block =np.zeros(shape=len(rds_demod_coeff))	#For Demodulation Resampler
	rds_rrc_block = np.zeros(shape=len((rrc_coeff)))			#For RRC convolution

	pll_block = np.array([0.0, 0.0, 1.0, 0.0, 1.0, 0.0])
	prevCarrier = np.zeros(shape=5121)
	prevCarrier = np.zeros(shape=5120)

	# audio buffer that stores all the audio blocks
	audio_data = np.array([]) # used to concatenate filtered blocks (audio data)

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

		# GET RDS CHANNEL
		channel_filt, channel_block = block_convolution(channel_coeff, dummy_fm, channel_block)

		# LOWER PATH OF RDS CARRIER RECOVERY: MATCH DELAY OF BPF (USED IN MIXER)
		channel_Delay, delay_block = myAllPass(channel_filt, delay_block)

		# UPPER PATH OF RDS CARRIER RECOVERY:
		carrier_Input = sq_nonlinearity(channel_filt)
		carrier_filt, carrier_block = block_convolution(carrier_coeff, carrier_Input, carrier_block)
		recoveredRDS, pll_block = fmPll(carrier_filt, 114e3, rf_Fs/rf_decim, 0.5, 0.0, 0.005, pll_block)

        # RDS DEMODULATION
		mixedAudio = mixer(recoveredRDS, channel_Delay)
		demod_filt, rds_demod_block = rs_block_convolution(rds_demod_coeff, mixedAudio, rds_demod_block, demod_decim, demod_exp)
		demod_filt, rds_rrc_block = block_convolution(rrc_coeff, demod_filt, rds_rrc_block)
		samples, sampling_intervals = CDR(demod_filt, sps)

		#Generate Plots of Monopath
		if block_count >= 3 and block_count < 6:
			print(sampling_intervals)
			n = 100
			x1 = range(n)
			x2 = range(n,n+n)
			fig2, axs = plt.subplots(3)
			fig2.suptitle('State saving checking')
			axs[0].plot(range(len(demod_filt)), demod_filt, c='blue')
			#axs[0].plot(x1, prev1[(len(prev1)-n):], c='orange')
			axs[0].plot(sampling_intervals, np.zeros(shape = len(sampling_intervals)), marker="x", c='orange', markersize=15)
			axs[0].plot(samples, np.zeros(shape = len(samples)), marker="o")
			axs[0].set_title('', fontstyle='italic',fontsize='medium')
			axs[0].axhline(y = 0, color = 'r', linestyle = '-')

			axs[1].plot(x2, carrier_filt[:n], c='blue')
			axs[1].plot(x1, prev2[(len(prev2)-n):], c='orange')
			axs[1].set_title('Mixed Audio', fontstyle='italic',fontsize='medium')
			axs[1].axhline(y = 0, color = 'r', linestyle = '-')

			axs[2].plot(x2, recoveredRDS[:n], c='blue')
			axs[2].plot(x1, prevPLL[(len(prevPLL)) - n:], c='orange')
			axs[2].set_title('PLL', fontstyle='italic',fontsize='medium')
			axs[2].axhline(y = 0, color = 'r', linestyle = '-')

			plt.show()
			quit()

		prev1 = demod_filt
		prev2 = carrier_filt
		prevPLL = recoveredRDS
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
