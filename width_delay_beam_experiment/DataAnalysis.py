# -*- coding: utf-8 -*-
"""
Created on Thu Nov 05 12:48:19 2015
Last updated: Sat Nov 07 00:12:15 2015

@author: Kevin Multani けヴぃん　むるたに

DataAnalysis.py contains the classes:

1) Waveform
2) Analyzer

The idea behind this code is to be able to create a one-click-result
batch analysis script from RAW data to the resulting plots.

The data is collected from the FFAG (Fixed-Field Alternating Gradient)
accelerator based in Kumatori, Japan as part of the Kyoto
University Research Reactor Institute.

Data is read from an oscilloscope and is analyzed as part of beam studies.
See corresponding report.
"""

import numpy as np
from scipy.fftpack import fft
from scipy.fftpack import fftfreq
from scipy.integrate import cumtrapz
import lmfit

class Waveform(object):

    def __init__(self, true_freq = 1.568e6):
        '''
        Waveform contains the raw data and 'cleaned' data. It
        provides methods to load and clean the data. The instance attributes
        make it easy to access these datas (no data is lost, there are
        variables containing both the raw and cleaned data).

        instance attributes:

        - filename(string): name of data file
        - data_path(string): path to where all the data is stored
        - true_freq(string): the frequency of the beam when the experiment was conducted
        - LINAC_delay(string): the delay of the trigger of the LINAC
        - pulse_width(string): the pulse width of the beam (duty cycle)

        - raw_V(np.array): the raw voltage data
        - raw_t(np.array): the raw time data
        - raw_numturns(np.array): time data scaled to show the turn number

        - clean_V(np.array): the cleaned voltage data
        - clean_t(np.array): the cleaned time data
        - clean_numturns(np.array): cleaned time data scaled to show turn number
        '''

        #### NAMING PARAMETERS ####
        self.filename = None # attribute set in load()
        self.data_path = None # attribute set in load()

        #### EXPERIMENT PARAMETERS ####
        self.true_freq = true_freq
        self.LINAC_delay = None # attribute set in load()
        self.pulse_width = None # attribute set in load()

        #### RAW DATA ####
        self.raw_V = None # attribute set in load()
        self.raw_t = None # attribute set in load()
        self.raw_numturns = None # attribute set in load()

        #### CLEANED DATA ####
        self.clean_V = None # attribute set in clean()
        self.clean_t = None # attribute set in clean()
        self.clean_numturns = None # attribute set in clean()
        return

    #### Methods ####
    def load(self, path_to_file, filename):
        '''
        load's purpose is to load in the raw data.

        assumes:
        - valid path_to_file and filenames

        - that the data file is comma separated and contains time and voltage
          data in columns 3 and 4, respectively.

        arguments:
        - path_to_file(string): the path to the folder containing data

        - filename: the filename of the data

        returns:
        - nothing: sets the filename, data_path, and raw data instance attributes
        '''

        # Set filename and data_path instance attribute
        self.filename = filename
        self.data_path = path_to_file

        # Read raw data and store
        mat = np.genfromtxt(self.data_path + self.filename, delimiter=',',
                            usecols=(3,4))

        # Set the raw data instance attributes
        self.raw_t, self.raw_V = mat[:,0]-mat[0,0], mat[:,1]
        self.raw_numturns = self.true_freq*self.raw_t

        # Set the experiment parameter instance attributes
        nocsv, empty = self.filename.split('.csv')
        self.LINAC_delay, self.pulse_width = (float(nocsv.split('_')[2]),
                                              float(nocsv.split('_')[-1]))
        return

    def clean(self, trim_start=10000, trim_stop=50000,
              noise_filname='rf_noise.csv'):
        '''
        clean has two purposes: trim the data to specify the signal region,
        and remove the rf noise data from the raw voltage data.

        If None or 'None' is passed in for the noise_filename, the data is
        only trimmed.

        assumes:
        that the load() method has been already called and there is raw
        data to be cleaned.

        arguments:
        - trim_start(int): the start index (excluding indices strictly less)
        - trim_stop(int): the stop index (excluding indices greater than or equal)
        - noise_filename(string): filename containing the noise data

        returns:
        nothing -- only sets the clean data instance attributes
        '''

        if noise_filname is None or noise_filname == 'None':

            # Set the cleaned voltage data attribute;
            # exclude indices < trim_start and indices >= trim_stop
            self.clean_V = self.raw_V[trim_start:trim_stop]

        else:

            # Read noise data and store
            noise_dat = np.genfromtxt(self.data_path + noise_filname,
                                      delimiter=',', usecols=(4,))

            # Subtract noise and raw voltage data
            temp = self.raw_V - noise_dat

            # Set the cleaned and noiseless voltage data attribute;
            # exclude indices < trim_start and indices >= trim_stop
            self.clean_V = temp[trim_start:trim_stop]

        # Set the cleaned time and num_turns data instance attributes;
        # exclude indices < trim_start and indices >= trim_stop
        self.clean_t = self.raw_t[trim_start:trim_stop]
        self.clean_numturns = self.raw_numturns[trim_start:trim_stop]
        return

class Analyzer(object):

    def __init__(self, waveform):
        '''
        Analyzer takes a Waveform object as input. The Anazlyer class provides
        a means to apply the necessary analysis (as determined by me). All of
        the class attributes are meant to store data processing at each step.
        The structure of this class is such that calculate_envelope
        should be called before calculate_bunch_area_diff -- as instance
        attributes are set in these methods (and thus creates dependencies).

        instance attributes:

        - waveform(Waveform): a waveform passed in when an Analyzer object is initialized

        - envelope_V(np.array): points which trace out the positive (or negative)
                      envelope passed in via waveform.
        - envelope_t(np.array): the time-axis corresponding with envelope_V

        - ideal_V(np.array): this is the cleaned voltage subtracted by the envelope
        - ideal_t(np.array): the time-corresponding with ideal_V
        - ideal_numturns(np.array): the number-of-turns-axis corresponding with envelope_V

        - int_bunch_diff(np.array): this is the difference of the integrated ideal_V per
                          turn. Corresponding to the integral;
                          int_i^{i+1} V_ideal dn (for i=[start,..,max num of turns])
        - int_bunch_diff_x(np.array): the x-axis corresponding to int_bunch_diff

        - fourier_V(np.array): the fft of V_clean
        - fourier_f(np.array): the frequency-axis of fourier_V
        - range_mask(np.array): a mask where when you call fourier_f[range_mask]
                      it gives you the interval 0.5*f_0 <= f <= 1.5*f_0
                      where f_0 is the revolution frequency.

        - alpha(float): the damping factor corresponding to beam loss
        - P(float): the parameter corresponding to the complex power i.e the frequency
             spectrum behaviour around 0.5*f_0 <= f <= 1.5*f_0 (f_0 is the
             revolution frequency).
        - measured_revfreq(float): the revolution frequency as per the peak in the fft
        - max_turns(int): the maximum number of turns corresponding with the
                          waveform passed into this Analyzer
        - start_turn_num(int): the beginning turn number of the waveform passed in
                               (since it has been cleaned up and points have been
                               excluded, start_turn_num is almost never 0).

        - envelope_V_predict(np.array): the fitted model of envelope_V
        - envelope_t_predict(np.array): the time-axis corresponding with envelope_V_predict
        - envelope_numturns_predict(np.array): the turn-number axis coressponding to
                                               envelope_V_predict

        - int_bunch_diff_predict(np.array): the fitted model of int_bunch_diff
        - int_bunch_diff_x_predict(np.array): the x-axis corresponding to
                                    int_bunch_diff_predict (it defaults to
                                    turn-number, but can be translated into
                                    time easily.)
        '''

        #### INITIALIZE DATA ####
        self.waveform = waveform # set via __init__

        #### ANALYZED DATA ####

        # Envelope data
        self.envelope_V = None # set in calculate_envelope()
        self.envelope_t = None # set in calculate_envelope()

        # 'Ideal' is the clean voltage data subtracted by envelope
        self.ideal_V = None # set in calculate_envelope()
        self.ideal_t = None # set in calculate_envelope()
        self.ideal_numturns = None # set in calculate_envelope()

        # Array that contains the integrated V difference per turn
        self.int_bunch_diff = None # set in calculate_bunchArea_diff()
        self.int_bunch_diff_x = None # set in calculate_bunchArea_diff()

        # FFT applied to clean voltage data
        self.fourier_V = None # set in calculate_P()
        self.fourier_f = None # set in calculate_P()
        self.fourier_range_mask = None # set in calculate_P()

        # Features from the data
        self.alpha = None # set in calculate_bunchArea_diff()
        self.P = None # set in calculate_P()
        self.measured_revfreq = None # set in calculate_P()

        # The stop and start of this data set in turn number --
        # can convert into time in seconds by dividing either
        # the measured frequency or the 'true' frequency attribute
        # of the waveform.
        self.max_turns = None # set in calculate_bunchArea_diff()
        self.start_turn_num = None # set in calculate_bunchArea_diff()

        #### FIT DATA ####

        # y and x values of the fitted envelope function
        self.envelope_V_predict = None # set in calculate_envelope()
        self.envelope_t_predict = None # set in calculate_envelope()
        self.envelope_numturns_predict = None # set in calculate_envelope()

        # y and x values of the integrated V difference per turn function
        self.int_bunch_diff_predict = None # set in calculate_bunchArea_diff()
        self.int_bunch_diff_x_predict = None # set in calculate_bunchArea_diff()
        return

    #### Methods #####
    def calculate_envelope(self, neg = False):
        '''
        Calculate_envelope: calculates the positive envelop of the given
        waveform by default (neg = True, for negative envelope) and fits
        the envelop to a model function.

        assumes:
        - the waveform contains clean data -- i.e. waveform.clean(..) has
          been called.

        arguments:
        - neg(bool): a boolean flag which determines if the positive envelope or
                     negative envelope (default False -- positive envelope)
                     is needed.

        returns:
        - nothing: sets the envelope, fitted envelope, and ideal
                   instance attributes.
        '''
        ## This is a hack to get the envelope, since the data is
        ## high frequency -- I keep the lookahead parameter very low
        ## and the delta parameter 0, so that it finds all the positive
        ## peaks of individual pulse signals.
        ##
        ## Then I run the algorithm again, on the result of the first run
        ## the net effect is a smoothing out of the maxima and we end up
        ## with the envelope.
        ##
        ## This sometimes has a bit of noise in the final result -- but
        ## works most of the time.
        tempx, tempy = peak_detect(self.waveform.clean_t, self.waveform.clean_V)
        self.envelope_t, self.envelope_V = peak_detect(tempx, tempy)

        # Exlcude a certain number of points from the beginning
        # of the data
        fitStart = 16

        # Fit the envelope
        mod, params = envelope_fit(self.envelope_t[fitStart:],
                          self.envelope_V[fitStart:],
                          verbose=False)

        # Set the predicted envelope instance attribute
        self.envelope_V_predict = mod.eval(params, t=self.waveform.clean_t)

        # Set the ideal_V instance attribute
        # ideal_V is what the signal should look like if the ground
        # voltage wasn't increasing.
        self.ideal_V = self.waveform.clean_V - self.envelope_V_predict

        # Re-set the predicted envelope data
        # instances to match len(self.envelope_V)
        self.envelope_V_predict = mod.eval(params, t=self.envelope_t)
        self.envelope_t_predict = self.envelope_t
        self.envelope_numturns_predict = self.envelope_t_predict*self.waveform.true_freq

        # Exclude positive values (since the amplifier at KURRI is an
        # inverting amplifier, the ideal signal is purely negative).
        filtered = self.ideal_V <= 0

        # Re-set the ideal values
        self.ideal_V, self.ideal_t, self.ideal_numturns = (
            self.ideal_V[filtered], self.waveform.clean_t[filtered],
            self.waveform.clean_numturns[filtered])
        return


    def calculate_bunchArea_diff(self):
        '''
        Calculates delta(n_i) = q_area(n_j) - q_area(n_i)
        for j = i + 1,
        where q_area(n) = int_0^(n) [V_clean(m) - V_envelope(m)] dm.

        delta(n_i) represents the change in integrated charge per consective
        turn / period.

        assumes:
        - self.waveform contains proper cleaned voltage and n data

        - self.envelope_V has been fitted and self.envelope_V_predict is
          calculated with the same t-values as self.waveform.cleaned_t

        arguments:
        - none

        returns:
        - nothing: method sets the int_bunch_diff, int_bunch_diff_predicted
                   (and corresponding x-axes), max_turns, start_turn_num, and
                   alpha instance attributes.
        '''

        # Initialize an empty numpy array
        self.int_bunch_diff = np.asarray([])

        # Determine the number of maximum number of turns in this set of data
        self.max_turns = int(max(self.ideal_numturns))

        # Determine the starting turn number (since we are
        # working with cleaned data)
        self.start_turn_num = int(min(self.ideal_numturns))


        # Check if ideal_V has values
        if self.ideal_V is None:
            print "You must call the calculate_envelope function."
            raise ValueError

        # For i < num_turns_array and num_turns_array <= i+1
        # calculate the integral
        for i in range(self.start_turn_num, self.max_turns+1):

            # Finds the indices that obey:
            # i < num_turns_array and num_turns_array <= i+1
            totalmask = np.logical_and(np.less(i, self.ideal_numturns),
                                       np.less_equal(self.ideal_numturns, i+1))

            # Calculate int_i^{i+1} integrand(n) dn
            Vinted = cumtrapz(self.ideal_V[totalmask],
                              self.ideal_numturns[totalmask],
                              initial=0)[-1]

            # Append the result into instance attribute
            self.int_bunch_diff = np.append(self.int_bunch_diff, Vinted)

        # The x-axis corresponding to self.int_bunch_diff
        self.int_bunch_diff_x = range(self.start_turn_num, self.max_turns+1)

        # Set the x-axis for int_bunch_predict (instance attribute)
        self.int_bunch_diff_x_predict = self.int_bunch_diff_x[:-1]

        # Create the variable to be passed as the x-axis to the fitting
        # method. self.int_bunch_diff's last point is erroneous, so I
        # exclude it, which explains the [:-1]
        n = np.asarray(range(len(self.int_bunch_diff)))[:-1]

        # Create a variable to be passed into the fitting function
        toFitData = self.int_bunch_diff[:-1]
        self.int_bunch_diff_predict, params = int_bunch_fit(n,
                                                            toFitData,
                                                            verbose=False)
        # Set the alpha instance attribute
        self.alpha = params['alpha'].value
        return


    def calculate_P(self):
        '''
        Calculates the P parameter, which is simply an integral around
        the fundamental revolution frequency given by the fft.
        The integral is on the interval: 0.5*f_0 <= f <= 1.5*f_0, where
        f_0 is the revolution frequency.

        assumes:
        - self.waveform contains proper voltage and time data.

        arguments:
        none

        returns:
        - nothing: sets the fourier data and P instance attributes.
        '''

        # Calculate the FFT of cleaned volatge and time data
        # and set the fourier instance attributes
        self.fourier_f, self.fourier_V = apply_fft(self.waveform.clean_t,
                                                   self.waveform.clean_V)


        # Create a boolean mask of frequencies 0.5*f_0 <= f <= 1.5*f_0
        self.range_mask = np.logical_and(
        np.less_equal(self.waveform.true_freq*0.5, self.fourier_f),
        np.less_equal(self.fourier_f, self.waveform.true_freq*1.5)
        )

        # Find the index of the peak at the revolution frequency
        ind_max = self.fourier_V[self.range_mask].argmax()

        # Get the value of the revolution frequency and set the instance
        # instance attribute
        self.measured_revfreq = self.fourier_f[self.range_mask][ind_max]

        # Integrate the spectrum in the interval 0.5*f_0 <= f <= 1.5*f_0
        # and set the P instance attribute
        self.P = cumtrapz(self.fourier_V[self.range_mask],
                          self.fourier_f[self.range_mask])[-1]
        return


#### Helper Methods #####
def step_responseRLC_overdamped(t, a_0, a_1, a_2, alpha, w_d):
    '''
    Callable function used as a model for fitting. This function
    is the overdamped step response for an RLC second order circuit.

    assumes:
    nothing

    arguments:
    - t(np.array): the independent variable
    - a_0(float): the dc component of response
    - a_1(float): related to f(0)
    - a_2(float): related to f'(0) and f(0)
    - alpha(float): the damping factor of the circuit
    - w_d(float): the damped frequency = sqrt(|alpha^2 - w_0^2|)

    returns:
    - the value at time t and parameter values a_0, a_1, a_2, w_d, t
    '''
    return (a_0 + a_1 * np.exp(-alpha * t) * np.cosh(w_d * t)
           + a_2 * np.exp(-alpha * t) * np.sinh(w_d * t))


def step_responseRLC_underdamped(t, a_0, a_1, a_2, alpha, w_d):
    '''
    Callable function used as a model for fitting. This function
    is the underdamped step response for an RLC second order circuit.

    assumes:
    nothing

    arguments:
    - t(np.array): the independent variable
    - a_0(float): the dc component of response
    - a_1(float): related to f(0)
    - a_2(float): related to f'(0) and f(0)
    - alpha(float): the damping factor of the circuit
    - w_d(float): the damped frequency = sqrt(|alpha^2 - w_0^2|)

    returns:
    - the value at time t and parameter values a_0, a_1, a_2, w_d, t
    '''
    return (a_0 + a_1 * np.exp(-alpha * t) * np.cos(w_d * t)
           + a_2 * np.exp(-alpha * t) * np.sin(w_d * t))


def apply_fft(time, voltage):
    '''
    Applies the FFT algorithm to the RAW voltage data.

    assumes:
    - len(time) == len(voltage)

    arguments:
    - time(np.array): a numpy array containing time sample points

    - voltage(np.array): a numpy array containing voltage sample points (at
      the corresponding times).

    returns:
    - xf(np.array): the frequency axis
    - yf(np.array): the magnitude of fft of voltage
    '''
    # Number of samplepoints
    N = len(voltage)

    # Sampling frequency
    Fs = np.diff(time).mean()

    # Apply FFT algorithm
    yf = fft(voltage)[0:N/2]

    # Normalize spectrum
    yf= np.abs(yf)

    # Create frequencies for x-axis
    xf = fftfreq(N, Fs)[:N/2]

    return xf, yf



def peak_detect(xarr, yarr, lookahead = 1, delta=0):
    '''
    Taken from :: sixtenbe (GitHub Gist)
    https://gist.github.com/sixtenbe

    Modified to be compatible with this file's structure.

    Converted from/based on a MATLAB script at:
    http://billauer.co.il/peakdet.html

    function for detecting local maximas and minmias in a signal.
    Discovers peaks by searching for values which are surrounded by lower
    or larger values for maximas and minimas respectively

    keyword arguments:
    y_arr -- A list containg the signal over which to find peaks
    x_arr -- (optional) A x-axis whose values correspond to the y_axis list
        and is used in the return to specify the postion of the peaks. If
        omitted an index of the y_axis is used. (default: None)
    lookahead -- (optional) distance to look ahead from a peak candidate to
        determine if it is the actual peak (default: 200)
        '(sample / period) / f' where '4 >= f >= 1.25' might be a good value
    delta -- (optional) this specifies a minimum difference between a peak and
        the following points, before a peak may be considered a peak. Useful
        to hinder the function from picking up false peaks towards to end of
        the signal. To work well delta should be set to delta >= RMSnoise * 5.
        (default: 0)
            delta function causes a 20% decrease in speed, when omitted
            Correctly used it can double the speed of the function

    return -- ym, xm : the maximum values and corresponding x values
    '''
    ## I don't need the minimum peaks feature, but the algorithm
    ## doesn't work as expected if I comment out anything that deals
    ## with the minimum. And I'm too lazy to fix it heh.
    maxPeaks = []
    minPeaks = []

    # Used to pop the first hit which almost always is false
    dump = []

    # Store data length for later use
    length = len(yarr)

    # Perform some checks
    if lookahead < 1:
        raise ValueError, "Lookahead must be '1' or above in value"
    if not (np.isscalar(delta) and delta >= 0):
        raise ValueError, "delta must be a positive number"

    # Maxima and minima candidates are temporarily stored in
    # mx and mn respectively
    mn, mx = np.Inf, -np.Inf

    # Only detect peak if there is 'lookahead' amount of points after it
    for index, (x, y) in enumerate(zip(xarr[:-lookahead],
                                        yarr[:-lookahead])):
        if y > mx:
            mx = y
            mxpos = x
        if y < mn:
            mn = y
            mnpos = x

        #### LOOK FOR MAX ####
        if y < mx-delta and mx != np.Inf:

            # Maxima peak candidate found
            # Look ahead in signal to ensure that this is a peak and not jitter
            if yarr[index:index+lookahead].max() < mx:
                maxPeaks.append([mxpos, mx])
                dump.append(True)

                # Set algorithm to only find minima now
                mx = np.Inf
                mn = np.Inf
                if index+lookahead >= length:
                    # End is within lookahead no more peaks can be found
                    break
                continue

        #### LOOK FOR MIN ####
        if y > mn+delta and mn != -np.Inf:

            # Minima peak candidate found
            # Look ahead in signal to ensure that this is a peak and not jitter
            if yarr[index:index+lookahead].min() > mn:
                minPeaks.append([mnpos, mn])
                dump.append(False)

                # Set algorithm to only find maxima now
                mn = -np.Inf
                mx = -np.Inf
                if index+lookahead >= length:

                    # End is within lookahead no more peaks can be found
                    break

    # Remove the false hit on the first value of the y_axis
    try:
        if dump[0]:
            maxPeaks.pop(0)
        else:
            minPeaks.pop(0)
        del dump
    except IndexError:
        # No peaks were found, do nothing
        pass

    ## Return only the maximum, that's all I need from this method
    xm, ym =  (np.asarray([p[0] for p in maxPeaks]),
               np.asarray([p[1] for p in maxPeaks]))
    return xm, ym


def envelope_fit(t, V, verbose = True):
    '''
    Uses the package lmfit to fit the given data: t and V. The model for this
    fit is an overdamped step response of a second order RLC circuit.

    assumes:
    - the data is suitable for the fitting model of overdamped step response
      of a second order RLC circuit.

    arguments:
    - t: the time-axis (or x-axis) of data
    - V: the voltage values

    returns:
    - fit_mod: the model object (see lmfit) corresponding to the model provided.
    - result.params: the Parameter (see lmfit) object corresponding to the
                     solution to the fitting algorithm.
    '''
    # Given initial parameters (set by eye after staring at the data for a while)
    init = [0.5, -2, 200, 1000, 50]

    # Set the model as a overdamped step response of RLC circuit
    fit_mod = lmfit.model.Model(step_responseRLC_overdamped)

    # Create the Parameter() object, pars is essentially a dictionary
    pars = fit_mod.make_params()

    # Set the parameters in the pars object to initial conditions.
    pars['a_0'].set(init[0])
    pars['a_1'].set(init[1])
    pars['a_2'].set(init[2])
    pars['alpha'].set(init[3])
    pars['w_d'].set(init[4])

    # Fit V and store resultant ModelResult object
    result = fit_mod.fit(V, pars, t=t)

    # verbose is True, print out the fit report (has stats info)
    if verbose:
        print(result.fit_report())

    return fit_mod, result.params

def int_bunch_fit(t, V, verbose = True):
    '''
    Uses the package lmfit to fit the given data: t and V. The model for this
    fit is an underdamped step response of a second order RLC circuit.

    assumes:
    - the data is suitable for the fitting model of overdamped step response
      of a second order RLC circuit.

    arguments:
    - t: the time-axis (or x-axis) of data
    - V: the voltage values

    returns:
    - result.best_fit: the y-values corresponding to the best fit (the best
                       parameters are in result.params)
    - result.params: the Parameter (see lmfit) object corresponding to the
                     solution to the fitting algorithm.
    '''

    # Given initial parameters (set by eye after staring at the data for a while)
    init = [-0.5, -0.5, 0.2, 0.055, 50]

    # Set the model as a underdamped step response of RLC circuit
    fit_mod = lmfit.model.Model(step_responseRLC_underdamped)

    # Create the Parameter() object, pars is essentially a dictionary
    pars = fit_mod.make_params()

    # Set the parameters in the pars object to initial conditions.
    pars['a_0'].set(init[0])
    pars['a_1'].set(init[1])
    pars['a_2'].set(init[2])
    pars['alpha'].set(init[3])
    pars['w_d'].set(init[4])

    # Fit V and store resultant ModelResult object
    result = fit_mod.fit(V, pars, t=t)

    # verbose is True, print out the fit report (has stats info)
    if verbose:
        print(result.fit_report())

    return result.best_fit, result.params
