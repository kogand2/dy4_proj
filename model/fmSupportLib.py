#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import numpy as np
import math, cmath

#
# you should use the demodulator based on arctan given below as a reference
#
# in order to implement your OWN FM demodulator without the arctan function,
# a very good and to-the-point description is given by Richard Lyons at:
#
# https://www.embedded.com/dsp-tricks-frequency-demodulation-algorithms/
#
# the demodulator boils down to implementing equation (13-117) from above, where
# the derivatives are nothing else but differences between consecutive samples
#
# needless to say, you should not jump directly to equation (13-117)
# rather try first to understand the entire thought process based on calculus
# identities, like derivative of the arctan function or derivatives of ratios
#

#
# use the four quadrant arctan function for phase detect between a pair of
# IQ samples; then unwrap the phase and take its derivative to demodulate
#

#======================================= Demodualtion for RL =======================================

def fmDemodArctan(I, Q, prev_phase = 0):
# the default prev_phase phase is assumed to be zero, however
# take note in block processing it must be explicitly controlled

	# empty vector to store the demodulated samples
	fm_demod = np.empty(len(I))

	# iterate through each of the I and Q pairs
	for k in range(len(I)):

		# use the atan2 function (four quadrant version) to detect angle between
		# the imaginary part (quadrature Q) and the real part (in-phase I)
		current_phase = math.atan2(Q[k], I[k])

		# we need to unwrap the angle obtained in radians through arctan2
		# to deal with the case when the change between consecutive angles
		# is greater than Pi radians (unwrap brings it back between -Pi to Pi)
		[prev_phase, current_phase] = np.unwrap([prev_phase, current_phase])

		# take the derivative of the phase
		fm_demod[k] = current_phase - prev_phase

		# save the state of the current phase
		# to compute the next derivative
		prev_phase = current_phase

	# return both the demodulated samples as well as the last phase
	# (the last phase is needed to enable continuity for block processing)
	return fm_demod, prev_phase

def myDemod(I, Q, dummy_state):
	fm_demod = np.empty(len(I)) #demodulated samples
	#insert prvious values into I and Q
	I_slice,Q_slice = np.split(dummy_state, 2)

	newI = np.concatenate((I_slice,I)) #Adds previous two values from the last block to array
	newQ = np.concatenate((Q_slice,Q)) #Adds previous two values from the last block to array

	for j in range(len(I)):
		#j + 2 is current value
		#j + 1 is the the previous value one index before the current values
		#j is the the previous value two indexs before the current values

		#Calculate the derivatives of I and Q then find the difference between them
		derQ = (newQ[j+1] - newQ[j])
		derI = (newI[j+1] - newI[j])
		fm_demod[j] = derQ*newI[j+1] - derI*newQ[j+1]
		#Scaling operation
		if ((newQ[j+1]**2 + newI[j+1]**2) != 0):
			fm_demod[j] = fm_demod[j]/(newQ[j+1]**2 + newI[j+1]**2)
		else:
			fm_demod[j] = 0
	dummy_state = np.array([I[j],Q[j]])
	#print(dummy_state)
	return fm_demod, dummy_state
# custom function for DFT that can be used by the PSD estimate

#======================================= filters =======================================
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

def myAllPass(input_block, state_block):

	output_block = np.concatenate((state_block, input_block[:-len(state_block)]))
	state_block = input_block[-len(state_block):]

	return output_block, state_block

def rrc(Fs, N_taps):

	"""
	Root raised cosine (RRC) filter

	Fs  		sampling rate at the output of the resampler in the RDS path
				sampling rate must be an integer multipler of 2375
				this integer multiple is the number of samples per symbol

	N_taps  	number of filter taps

	"""

	# duration for each symbol - should NOT be changed for RDS!
	T_symbol = 1/2375.0

	# roll-off factor (must be greater than 0 and smaller than 1)
	# for RDS a value in the range of 0.9 is a good trade-off between
	# the excess bandwidth and the size/duration of ripples in the time-domain
	beta = 0.90

	# the RRC inpulse response that will be computed in this function
	impulseResponseRRC = np.empty(N_taps)

	for k in range(N_taps):
		t = float((k-N_taps/2))/Fs
		# we ignore the 1/T_symbol scale factor
		if t == 0.0: impulseResponseRRC[k] = 1.0 + beta*((4/math.pi)-1)
		elif t == -T_symbol/(4*beta) or t == T_symbol/(4*beta):
			impulseResponseRRC[k] = (beta/np.sqrt(2))*(((1+2/math.pi)* \
					(math.sin(math.pi/(4*beta)))) + ((1-2/math.pi)*(math.cos(math.pi/(4*beta)))))
		else: impulseResponseRRC[k] = (math.sin(math.pi*t*(1-beta)/T_symbol) +  \
					4*beta*(t/T_symbol)*math.cos(math.pi*t*(1+beta)/T_symbol))/ \
					(math.pi*t*(1-(4*beta*t/T_symbol)*(4*beta*t/T_symbol))/T_symbol)

	# returns the RRC impulse response to be used by convolution
	return impulseResponseRRC
#======================================= Convolution and sampling =======================================

def block_convolution(h, xb, state): #Edits made from Feb 10, partitioned design
	yb = np.zeros(len(xb))
	stateLen = len(state)

	for n in range(len(yb)):
		for k in range(len(h)):
			if n-k >= 0:
				yb[n] += h[k] * xb[n-k]
			else:
				yb[n] += h[k] * state[stateLen + (n - k)]

	new_state = xb[len(xb) - len(h) + 1:]

	return yb, new_state

def ds_block_convolution(h, x, state, decim):
	y = np.zeros(len(x))
	down =  np.array([])
	stateLen = len(state)

	#count = 0;

	for n in range(0, len(x), decim):		#dominant partition
		y[n] = 0.0
		for k in range(len(h)):
			if (n-k) >= 0:
				y[n] += x[n-k]*h[k]
			else:
				y[n] += state[stateLen + n - k]*h[k]
			#count++
		down = np.append(down, y[n])

	new_state = x[len(x) - len(h) + 1:]

	return down, new_state

def rs_block_convolution(h, x, state, decim, exp):
	y = np.zeros(len(x)*exp)
	resample =  np.array([])
	#x_index = 0

	for n in range(0, len(x)*exp, decim):		#dominant partition
		y[n] = 0.0
		phase = n % exp
		x_index = int((n-phase)/exp)

		for k in range(phase,len(h),exp):
			if x_index >= 0:
				y[n] += x[x_index]*h[k]*exp
			else:
				y[n] += exp * h[k] * state[len(state) + x_index]

			x_index -= 1

		resample = np.append(resample, y[n])

	new_state = x[len(x) - len(h)//exp + 1:]

	return resample, new_state
#======================================= path specific =======================================

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
	ncoOutI = np.empty(len(pllIn)+1)
	ncoOutQ = np.empty(len(pllIn)+1)

	# initialize internal state
	integrator = state[0]
	phaseEst = state[1]
	feedbackI = state[2]
	feedbackQ = state[3]
	ncoOutI[0] = state[4]
	ncoOutQ[0] = state[5]
	trigOffset = state[6]
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
		ncoOutI[k+1] = math.cos(trigArg*ncoScale + phaseAdjust)
		ncoOutQ[k+1] = math.sin(trigArg*ncoScale + phaseAdjust)
	#print(trigOffset)
	new_state = [integrator, phaseEst, feedbackI, feedbackQ, ncoOutI[-1], ncoOutQ[-1], trigOffset]

	# for stereo only the in-phase NCO component should be returned
	# for block processing you should also return the state
	#print(ncoOut)
	return ncoOutI, ncoOutQ, new_state
	# for RDS add also the quadrature NCO component to the output

def mixer(recoveredStereo, channel_filt):
	mixedAudio = np.empty(len(channel_filt))
	for i in range(len(channel_filt)):
		mixedAudio[i] = 2 * recoveredStereo[i] * channel_filt[i]	#this would have both the +ve and -ve part of the cos combined, we need to keep the -ve part and filter it
																					#could prob automatically be done using the filter code from mono (???)
	return mixedAudio

def sq_nonlinearity(signalIn):
	signalLen = len(signalIn)
	output = np.zeros(shape = signalLen)
	for i in range(signalLen):
		output[i] = signalIn[i]*signalIn[i]

	return output #remember output will have a squared half amplitude and a positve offset of the same mag

def diff_decoding(manchester_values, initial, cdr_state):
	bits = []
	decoded_bits = []
	is_first_block = False

	if (initial == -1):
		is_first_block = True
		ind1 = 0
		ind2 = 0

		for i in range(0, len(manchester_values) - 1, 2):
			if manchester_values[i] == 1 and manchester_values[i+1] == 1: # HH (ignore)
				ind1 += 1
			elif manchester_values[i] == 0 and manchester_values[i+1] == 0: # LL (ignore)
				ind1 += 1

		for i in range(1, len(manchester_values) - 1, 2):
			if manchester_values[i] == 1 and manchester_values[i+1] == 1: # HH (ignore)
				ind2 += 1
			elif manchester_values[i] == 0 and manchester_values[i+1] == 0: # LL (ignore)
				ind2 += 1

		if (ind1 > ind2):
			initial = 1
		else:
			initial = 0

	# in the case that we need state-saving we use this index
	# which will NOT be returned and end up overwriting the initial index
	# from block 1
	start_index = initial
	if initial == 1 and not is_first_block:
		#print("sample_val[-1]: ",cdr_state[0])
		manchester_values.insert(0, cdr_state[0])
		start_index = 0

	#print("CHOSE : " +str(start_index))
	#print(manchester_values)
	for i in range(start_index, len(manchester_values) - 1, 2):
		if manchester_values[i] == 0 and manchester_values[i+1] == 1:	# LH = 0
			bits.append(0)
		elif manchester_values[i] == 1 and manchester_values[i+1] == 0: # HL = 1
			bits.append(1)
		# change the two H's into one H
		elif manchester_values[i] == 1 and manchester_values[i+1] == 1: # HH (ignore)
			#i -= 1
			#print("Consecutive HI")
			pass
		# change the two L's into one L
		elif manchester_values[i] == 0 and manchester_values[i+1] == 0: # LL (ignore)
			#i -= 1
			#print("Consecutive LO")
			pass


	if not is_first_block:
		bits.insert(0, cdr_state[2])

	else:
		decoded_bits.append(bits[0])

	for i in range(1, len(bits)):
		decoded_bits.append(bits[i] ^ bits[i-1])
		#print(i)


	#print("THIS IS MANCHESTER")
	#print(manchester_values)
	#print("Manchester Length: " + str(len(manchester_values)))

	#print("THIS IS CDR State")
	#print(cdr_state)
	#print("Manchester Length: " + str(len(cdr_state)))

	#print("THIS IS DECODED")
	#print(decoded_bits)
	#print("Decoded Length: " + str(len(decoded_bits)))
	#print("===================================================================")
	return decoded_bits, initial, bits[-1]

def CDR_state(signalIn, interval, initial):

	samples = []				# x vals (just for testing)
	sample_vals = []			# y vals (manchester_values)
	best_init = initial			# index sampling starts from

	# processing first block, find best initial point to start sampling
	if (initial == -1):
		max = 0.0
		curr = 0.0
		best_init = 0

		for initial in range(0, interval):
			# reset for next vals
			curr = 0.0
			for i in range(initial, len(signalIn), interval):
				curr += abs(signalIn[i])

			#print("For initial point = " + str(initial) + " : " + str(curr))
			if (max < curr):
				max = curr
				best_init = initial

	#print("Chose initial " + str(best_init))
	# processing rest of blocks, finding values of samples
	for i in range(best_init, len(signalIn), interval):
		samples.append(i)
		if (signalIn[i] > 0):	# get a HI
			sample_vals.append(1)
		else:					# get a LO
			sample_vals.append(0)

	#print("===================================================================")
	return samples, sample_vals, best_init

def frame_sync(decoded_bits, prev_decoded):
	start_point = 0 #relative to bit stream
	if len(prev_decoded) == 0:
		return None, start_point, decoded_bits
	else:
		#print("This is the obtained bit_stream")
		bit_stream = np.append(prev_decoded, decoded_bits)
		#print(bit_stream)
		#print(len(bit_stream))
	if len(bit_stream) >=  26:
		block = ["A","B","C","D"]
		for i in range(len(bit_stream) - 25):
			check = bit_stream[i:26+i]
			message = matrixMult(check, getParityCheck())
			for k in range(4):
				if np.array_equal(message, getSyndrome(k)):
					print("CORRECT SYNDROME OBTAINED")
					print(check)
					print("obtained: " + block[k] + " " + str(i))
					return block[k], i, prev_decoded


	return None, start_point, bit_stream[len(bit_stream) - 25:]

def app_layer(block_type, start_point, decoded_bits, prev_decoded, d_service, d_index):

	#print("This is the obtained bit_stream")
	bit_stream = np.append(prev_decoded, decoded_bits)
	#print("Bit stream, previous and decoded bits")
	#print(bit_stream)
	#print(len(bit_stream))
	#print(prev_decoded)
	#print(decoded_bits)
	#print(len(decoded_bits))
	#print(len(bit_stream))
	print("")
	if len(bit_stream) >=  26: #Repition 11.4 per second
		print("Block obtained")
		#print(bit_stream[start_point:26+start_point])
		if(block_type == "A"):
			#print("(====~~~~~~~~~~BLOCK A~~~~~~~~~~====)")
			PI_hex = ''
			for i in range(16):
				PI_hex+= str(bit_stream[start_point+i])

			#print("ORIGINAL PI IS ", PI_hex)

			PI_hex.replace('0x','')

			#print("The PI code in binary is: ", PI_hex)

			PI_hex = hex(int(PI_hex, 2))

			#print("The PI code is: ", PI_hex[2:])

			#print("Now entering block B")
			block_type = "B"

		elif(block_type == "B"):
			#Unsure if masking will be needed
			#print("(====~~~~~~~~~~BLOCK B~~~~~~~~~~====)")
			if np.array_equal(bit_stream[start_point:start_point+5],np.array([0,0,0,0,0])):
				#print("Group type is ", bit_stream[start_point:start_point+5])
				#print("Program Type ", bit_stream[start_point+6:start_point+11])
				d_index = int(str(bit_stream[start_point+14])+str(bit_stream[start_point+15]), 2)

			else:
				pass
				#print("ERROR: type is not 0A")

			#print("d_index is " + str(d_index))
			#print("Now entering block C")
			block_type = "C"

		elif(block_type == "C"):
			#print("(====~~~~~~~~~~BLOCK C~~~~~~~~~~====)")
			#print("--------------")
			#print("--------------")
			#print("Now entering block D")
			block_type = "D"

		elif(block_type == "D"):
			print("(====~~~~~~~~~~BLOCK D~~~~~~~~~~====)")
			d1 = d2 = ""
			arr = []
			for i in range(8):
				d1 += str(bit_stream[start_point+i])
				d2 += str(bit_stream[start_point+i+8])

			print(d1 + " " + d2)
			print("d_index is " + str(d_index))
			arr = bin_to_char(d1 + " " + d2)
			print("to be inserted is " + str(arr))
			if (len(d_service) < 8):

				if (d_index == 0 and len(d_service) == 0):
					d_service.append(arr[0])
					d_service.append(arr[1])

				elif (d_index == 1 and len(d_service) == 2):
					d_service.append(arr[0])
					d_service.append(arr[1])

				elif (d_index == 2 and len(d_service) == 4):
					d_service.append(arr[0])
					d_service.append(arr[1])

				elif (d_index == 3 and len(d_service) == 6):
					d_service.append(arr[0])
					d_service.append(arr[1])

			elif (d_index == 0):
				print("Program Service: " + "".join(d_service))
				d_service = []
				d_service.append(arr[0])
				d_service.append(arr[1])
			#print(d_service)
			#print("Now entering block A")
			block_type = "A"
		else:
			print("Unknown block type")

		return block_type, 0, bit_stream[start_point+26:], d_service, d_index

	else:
		start_point = 0 #relative to bit stream
		return block_type, start_point, bit_stream, d_service, d_index

def bin_to_char(a_binary_string):
	binary_values = a_binary_string.split()

	ascii_string = []
	for binary_value in binary_values:
	    an_integer = int(binary_value, 2)
	    ascii_character = chr(an_integer)
	    ascii_string.append(ascii_character)


	return ascii_string

def matrixMult(x, y):
	# This is not a generic matrix multiplication
	# This is build specificaly for bits and expected dimensions

	result = np.array([])
	value = 0
	# iterate through columns of y
	for i in range(len(y[0])):# iterate through columns of x and rows of y

		for k in range(len(x)):
			#print(str(x[k]) + " and " + str(y[i][k]))
			#print(k)
			if k == 0:
				value = (int(x[k]) and y[k][i])
			else:
				value = value ^ (int(x[k]) and y[k][i])
		result = np.append(result, value)
	#print(result)
	return result

def getParityCheck():
	parityCheck = np.array([[1,0,0,0,0,0,0,0,0,0],\
	[0,1,0,0,0,0,0,0,0,0],\
	[0,0,1,0,0,0,0,0,0,0],\
	[0,0,0,1,0,0,0,0,0,0],\
	[0,0,0,0,1,0,0,0,0,0],\
	[0,0,0,0,0,1,0,0,0,0],\
	[0,0,0,0,0,0,1,0,0,0],\
	[0,0,0,0,0,0,0,1,0,0],\
	[0,0,0,0,0,0,0,0,1,0],\
	[0,0,0,0,0,0,0,0,0,1],\
	[1,0,1,1,0,1,1,1,0,0],\
	[0,1,0,1,1,0,1,1,1,0],\
	[0,0,1,0,1,1,0,1,1,1],\
	[1,0,1,0,0,0,0,1,1,1],\
	[1,1,1,0,0,1,1,1,1,1],\
	[1,1,0,0,0,1,0,0,1,1],\
	[1,1,0,1,0,1,0,1,0,1],\
	[1,1,0,1,1,1,0,1,1,0],\
	[0,1,1,0,1,1,1,0,1,1],\
	[1,0,0,0,0,0,0,0,0,1],\
	[1,1,1,1,0,1,1,1,0,0],\
	[0,1,1,1,1,0,1,1,1,0],\
	[0,0,1,1,1,1,0,1,1,1],\
	[1,0,1,0,1,0,0,1,1,1],\
	[1,1,1,0,0,0,1,1,1,1],\
	[1,1,0,0,0,1,1,0,1,1]])
	return parityCheck

def getOffset(row):
	#Offset array
	offset = np.array([[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0],\
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,0,0],\
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,1,0,0,0],\
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,1,0,0]])

	return offset[row]

def getSyndrome(row):
	#expected Syndrome array
	expSyndrome = np.array([[1,1,1,1,0,1,1,0,0,0],\
	[1,1,1,1,0,1,0,1,0,0],\
	[1,0,0,1,0,1,1,1,0,0],\
	[1,0,0,1,0,1,1,0,0,0]])

	return expSyndrome[row]

#======================================= from lab 3 =======================================
def DFT(x):

	# number of samples
	N = len(x)

	# frequency bins
	Xf = np.zeros(N, dtype='complex')

	# iterate through all frequency bins/samples
	for m in range(N):
		for k in range(N):
			Xf[m] += x[k] * cmath.exp(1j * 2 * math.pi * ((-k) * m) / N)

	# return the vector that holds the frequency bins
	return Xf

# custom function to estimate PSD based on the Bartlett method
# this is less accurate than the Welch method from matplotlib
# however, as the visual inspections confirm, the estimate gives
# the user a "reasonably good" view of the power spectrum
def estimatePSD(samples, NFFT, Fs):

	# rename the NFFT argument (notation consistent with matplotlib.psd)
	# to freq_bins (i.e., frequency bins for which we compute the spectrum)
	freq_bins = NFFT
	# frequency increment (or resolution of the frequency bins)
	df = Fs/freq_bins

	# create the frequency vector to be used on the X axis
	# for plotting the PSD on the Y axis (only positive freq)
	freq = np.arange(0, Fs/2, df)

	# design the Hann window used to smoothen the discrete data in order
	# to reduce the spectral leakage after the Fourier transform
	hann = np.empty(freq_bins)
	for i in range(len(hann)):
		hann[i] = pow(math.sin(i*math.pi/freq_bins),2)
		#print(hann[i])
	# create an empty list where the PSD for each segment is computed
	psd_list = []

	# samples should be a multiple of frequency bins, so
	# the number of segments used for estimation is an integer
	# note: for this to work you must provide an argument for the
	# number of frequency bins not greater than the number of samples!
	no_segments = int(math.floor(len(samples)/float(freq_bins)))

	# iterate through all the segments
	for k in range(no_segments):

		# apply the hann window (using pointwise multiplication)
		# before computing the Fourier transform on a segment
		windowed_samples = samples[k*freq_bins:(k+1)*freq_bins] * hann
		#print(windowed_samples)
		# compute the Fourier transform using the built-in FFT from numpy
		Xf = np.fft.fft(windowed_samples, freq_bins)

		# note, you can check how MUCH slower is DFT vs FFT by replacing the
		# above function call with the one that is commented below
		#
		#Xf = DFT(windowed_samples)

		#
		# note: the slow impelementation of the Fourier transform is not as
		# critical when computing a static power spectra when troubleshooting
		#
		# note also: time permitting a custom FFT can be implemented


		# since input is real, we keep only the positive half of the spectrum
		# however, we will also add the signal energy of negative frequencies
		# to have a better a more accurate PSD estimate when plotting
		Xf = Xf[0:int(freq_bins/2)] # keep only positive freq bins
		#print(Xf)
		#print(freq_bins/2)
		#print(len(Xf))
		psd_seg = (1/(Fs*freq_bins/2)) * (abs(Xf)**2) # compute signal power
		psd_seg = 2*psd_seg # add the energy from the negative freq bins

		# translate to the decibel (dB) scale
		for i in range(len(psd_seg)):
			psd_seg[i] = 10*math.log10(psd_seg[i])


		# append to the list where PSD for each segment is stored
		# in sequential order (first segment, followed by the second one, ...)
		psd_list.extend(psd_seg)
		#print(k)
	# compute the estimate to be returned by the function through averaging
	psd_est = np.zeros(int(freq_bins/2))

	# iterate through all the frequency bins (positive freq only)
	# from all segments and average them (one bin at a time ...)
	for k in range(int(freq_bins/2)):
		# iterate through all the segments
		for l in range(no_segments):
			psd_est[k] += psd_list[k + l*int(freq_bins/2)]
		# compute the estimate for each bin
		psd_est[k] = psd_est[k] / no_segments

	# the frequency vector and PSD estimate
	return freq, psd_est

# custom function to format the plotting of the PSD
def fmPlotPSD(ax, samples, Fs, height, title):

	x_major_interval = (Fs/12)		# adjust grid lines as needed
	x_minor_interval = (Fs/12)/4
	y_major_interval = 20
	x_epsilon = 1e-3
	x_max = x_epsilon + Fs/2		# adjust x/y range as needed
	x_min = 0
	y_max = 10
	y_min = y_max-100*height
	ax.psd(samples, NFFT=512, Fs=Fs)
	#
	# below is the custom PSD estimate, which is based on the Bartlett method
	# it less accurate than the PSD from matplotlib, however it is sufficient
	# to help us visualize the power spectra on the acquired/filtered data
	#
	# freq, my_psd = estimatePSD(samples, NFFT=512, Fs=Fs)
	# ax.plot(freq, my_psd)
	#
	ax.set_xlim([x_min, x_max])
	ax.set_ylim([y_min, y_max])
	ax.set_xticks(np.arange(x_min, x_max, x_major_interval))
	ax.set_xticks(np.arange(x_min, x_max, x_minor_interval), minor=True)
	ax.set_yticks(np.arange(y_min, y_max, y_major_interval))
	ax.grid(which='major', alpha=0.75)
	ax.grid(which='minor', alpha=0.25)
	ax.set_xlabel('Frequency (kHz)')
	ax.set_ylabel('PSD (db/Hz)')
	ax.set_title(title)

if __name__ == "__main__":

	# do nothing when this module is launched on its own

	pass
