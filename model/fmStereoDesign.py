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
from fmSupportLib import fmDemodArctan, fmPlotPSD, myDemod
#from fmMonoBlock import *
# for take-home add your functions
def myLowPass(fc, fs, Ntaps):
	normCutoff = fc/(fs/2)
	h = np.zeros(Ntaps)
	for i in range(Ntaps):
		if (i==(Ntaps-1)/2):
			h[i] = normCutoff
		else:
			h[i] = normCutoff*(math.sin(np.pi*normCutoff*(i-(Ntaps-1)/2)))/(np.pi*normCutoff*(i-(Ntaps-1)/2))
		h[i] = h[i]*(math.sin(i*np.pi/Ntaps))**2

	return  h

def myBandPass(fb, fe, fs, Ntaps):
	normCenter = ((fe + fb)/2)/(fs/2)
	normPass = (fe - fb)/(fs/2)
	h = np.zeros(Ntaps)

	for i in range(Ntaps):
		if (i==(Ntaps-1)/2):
			h[i] = normPass # avoid division by zero in sinc for the center tap when Ntaps is odd
		else:
			h[i] = normPass*((math.sin(np.pi*(normPass/2)*(i-(Ntaps-1)/2)))/(np.pi*(normPass/2)*(i-(Ntaps-1)/2)))
		h[i] = h[i]*(math.cos(i*np.pi*normCenter))
		h[i] = h[i]*((math.sin(i*np.pi/Ntaps))**2)

	return  h

def block_convolution(h, xb, state):
	yb = np.zeros(len(xb))
	stateLen = len(state)

	for n in range(len(yb)):
		for k in range(len(h)):
			if n-k >= 0:
				yb[n] += h[k] * xb[n-k]
			else:
				yb[n] += h[k] * state[stateLen + (n - k)]
	#print(stateLen)
	#print(len(xb))
	#print(len(h))
	new_state = xb[len(xb) - len(h) + 1:]

	return yb, new_state

def fmPll(pllIn, freq, Fs, ncoScale = 2.0, phaseAdjust = 0.0, normBandwidth = 0.01, state = []):

	"""
	pllIn 	 		array of floats
					input signal to the PLL (assume known frequency)

	freq 			float
					reference frequency to which the PLL locks

	Fs  			float
					sampling rate for the input/output signals

	ncoScale		float
					frequency scale factor for the NCO output

	phaseAdjust		float
					phase adjust to be added to the NCO output only

	normBandwidth	float
					normalized bandwidth for the loop filter
					(relative to the sampling rate)

	state 			to be added

	"""

	# scale factors for proportional/integrator terms
	# these scale factors were derived assuming the following:
	# damping factor of 0.707 (1 over square root of 2)
	# there is no oscillator gain and no phase detector gain
	Cp = 2.666
	Ci = 3.555

	# gain for the proportional term
	Kp = (normBandwidth)*Cp
	# gain for the integrator term
	Ki = (normBandwidth*normBandwidth)*Ci

	# output array for the NCO
	ncoOut = np.empty(len(pllIn)+1)

	# initialize internal state
	integrator = state[0]
	phaseEst = state[1]
	feedbackI = state[2]
	feedbackQ = state[3]
	ncoOut[0] = state[4]
	trigOffset = state[5]
	# note: state saving will be needed for block processing


	for k in range(len(pllIn)):

		# phase detector
		errorI = pllIn[k] * (+feedbackI)  # complex conjugate of the
		errorQ = pllIn[k] * (-feedbackQ)  # feedback complex exponential

		# four-quadrant arctangent discriminator for phase error detection
		errorD = math.atan2(errorQ, errorI)

		# loop filter
		integrator = integrator + Ki*errorD

		# update phase estimate
		phaseEst = phaseEst + Kp*errorD + integrator


		# internal oscillator
		trigOffset += 1.0
		trigArg = (2*np.pi*(freq/Fs)*(trigOffset) + phaseEst)
		feedbackI = math.cos(trigArg)
		feedbackQ = math.sin(trigArg)
		ncoOut[k+1] = math.cos(trigArg*ncoScale + phaseAdjust)
	#print(trigOffset)
	new_state = [integrator, phaseEst, feedbackI, feedbackQ, ncoOut[-1], trigOffset]

	# for stereo only the in-phase NCO component should be returned
	# for block processing you should also return the state
	#print(ncoOut)
	return ncoOut, new_state
	# for RDS add also the quadrature NCO component to the output

def mixer(recoveredStereo, channel_filt): #Possible issues over here
	mixedAudio = np.empty(len(channel_filt))
	for i in range(len(channel_filt)):
		mixedAudio[i] = 2 * recoveredStereo[i] * channel_filt[i]	#this would have both the +ve and -ve part of the cos combined, we need to keep the -ve part and filter it
																					#could prob automatically be done using the filter code from mono (???)
	return mixedAudio

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
	filt_block = np.zeros(shape=len(carrier_coeff))		#For Stereo Processing
	carrier_block = np.zeros(shape=len(carrier_coeff))	#For Stereo Carrier Recovery
	channel_block = np.zeros(shape=len(channel_coeff))	#For Stereo Channel Extraction
	mono_block = np.zeros(shape=len(audio_coeff))		#For Monopath

	left_block = np.zeros(shape=1024)		#For Monopath
	right_block = np.zeros(shape=1024)		#For Monopath

	pll_block = np.array([0.0, 0.0, 1.0, 0.0, 1.0, 0.0])
	# add state as needed for the mono channel filter

	prevStereo = np.zeros(shape=5121)
	prevCarrier = np.zeros(shape=5120)

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
		dummy_fm, dummy_state = myDemod(i_ds, q_ds, dummy_state)

    #RF front end stops here==================================================

		# extract the stereo audio data

        #Stereo Carrier Recovery: Bandpass -> PLL -> Numerically Controlled Oscillator
		carrier_filt, carrier_block = block_convolution(carrier_coeff, dummy_fm, carrier_block)
		recoveredStereo, pll_block = fmPll(carrier_filt, 19e3, rf_Fs/rf_decim, 2.0, 0.0, 0.01, pll_block)

        #Stereo Channel Extraction: Bandpass
		channel_filt,channel_block = block_convolution(channel_coeff, dummy_fm, channel_block)

        #Stereo Processing: Mixer -> Digital filtering (Lowpass -> down sample) -> Stereo Combiner
		mixedAudio = mixer(recoveredStereo, channel_filt)

		stereo_filt,filt_block = block_convolution(audio_coeff, mixedAudio, filt_block)

		stereo_block = stereo_filt[::audio_decim]

		# extract the mono audio data
		audio_filt,mono_block = block_convolution(audio_coeff, dummy_fm, mono_block)

		audio_block = audio_filt[::audio_decim]

		for i in range(len(audio_block)):
			left_block[i] = audio_block[i] + stereo_block[i]

		for i in range(len(audio_block)):
			right_block[i] = audio_block[i] - stereo_block[i]

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

		if block_count >= 10 and block_count < 12:

			# plot PSD of selected block after FM demodulation
			ax0.clear()
			fmPlotPSD(ax0, dummy_fm, (rf_Fs/rf_decim)/1e3, subfig_height[0], \
					'Demodulated FM (block ' + str(block_count) + ')')
			# output binary file name (where samples are written from Python)
			fm_demod_fname = "../data/fm_demod_" + str(block_count) + ".bin"
			# create binary file where each sample is a 32-bit float
			dummy_fm.astype('float32').tofile(fm_demod_fname)

			# plot PSD of selected block after extracting mono audio
			# ... change as needed
			fmPlotPSD(ax1, channel_filt, ((rf_Fs/rf_decim)/1e3), subfig_height[1], 'channel_filt')

			fmPlotPSD(ax2, right_data, ((rf_Fs/rf_decim)/1e3), subfig_height[2], 'Right Channel')
			# plot PSD of selected block after downsampling mono audio
			# ... change as needed
			fmPlotPSD(ax3, left_data, audio_Fs/1e3, subfig_height[3], 'Left Channel')
			# save figure to file
			fig.savefig("../data/fmMonoBlock" + str(block_count) + ".png")

		block_count += 1

	print('Finished processing all the blocks from the recorded I/Q samples')

	complete_data = np.vstack((left_data, right_data)).T

	# write audio data to file
	out_fname = "../data/fmStereoBlock.wav"
	wavfile.write(out_fname, int(audio_Fs), np.int16((complete_data/2)*32767))
	print("Written audio samples to \"" + out_fname + "\" in signed 16-bit format")

	# uncomment assuming you wish to show some plots
	# plt.show()
