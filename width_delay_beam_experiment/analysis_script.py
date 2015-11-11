# -*- coding: utf-8 -*-
"""
Created on Sat Nov 07 14:11:38 2015

@author: Kevin Multani けヴぃん　むるたに

anaylsis_script.py performs analysis on all the data and creates the necessary
plots. All plots that are in the report are generated in this file. See
DataAnalysis.py for all the details of the analysis.
"""

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import DataAnalysis
import os

# Set all figures to render in LaTeX
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Path where the data is stored
data_path = u'C:\\Users\\Kevin\\Dropbox\\mylife\\coop\\Japan(offered,acccepted)\\work\\width_delay_beam_experiment\\data_26-10-2015\\'

# The dictionary where color hex values are the labels and keys are the type
# of plot
color_dict = {'data': u'#6699FF', 'fit': u'#00CC00', u'residuals': u'#FF6666'}

# The filename where noise data is contained
filename_noise = u'rf_noise.csv'

# The filename of the data that serves as the sample plots
filename_sample = u'inu_delay_1570_width_19.9.csv'

def plot_noise(path_to_data, filename, save_path='.//'):
    '''
    Plots the noise data with pretty formatting.

    assumes:
    - The data is a csv file where the time and voltage data are in columns
      3 and 4, respectively.

    - save_path, if given, is a valid directory.

    arguments:
    - path_to_data(string): The absolute or relative path to the
                            folder containing the noise data file.

    - filename(string): the filename of noise data file

    - save_path(string): the path where the plot gets saved (default to current
                         working directory)

    returns:
    - nothing: simply plots and saves the plot in the directory given by
               save_path.
    '''
    # Read data, the data is in columns 3 and 4 in a comma separated file
    dat = np.genfromtxt(path_to_data+filename, delimiter=',', usecols=(3,4))

    # Assign x and y values from data
    x, y = dat[:,0], dat[:,1]

    # Initialize figure and axis objects
    fig = plt.figure(figsize = (10,8))
    ax = plt.gca()

    # Plot the the data
    ax.plot(x, y, c=color_dict['data'], label=filename.replace('_', '\_'), alpha = 0.7)

    # Begin formatting plot

    # Change xtick and ytick frequency
    for label in ax.get_xticklabels()[1::2]:
        label.set_visible(False)
    for label in ax.get_yticklabels()[1::2]:
        label.set_visible(False)

    # Adding horizontal grid lines
        ax.yaxis.grid(True)

    # Remove axis spines
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(True)
    ax.spines["left"].set_visible(False)

    # Labels
    plt.xlabel('Time $[s]$', fontsize=15)
    plt.ylabel('Voltage $[V]$', fontsize=15)
    plt.legend(loc='best', fancybox=True, fontsize=15)

    # Tight layout
    plt.tight_layout()

    # Save the plot
    plt.savefig(save_path+filename.split('.')[0] + '.pdf')
    return

def plot_sample_waveform(analyzed_waveform, save_path='.//'):
    '''
    Plots and saves a total of four figures. The purpose of this function
    is to display the analysis process, step by step of a sample dataset.

    Plot 1): Plots the raw data and cleaned data
    Plot 2): Plots the cleaned data, the envelope of the cleaned data,
             the fit for the envelope and the residuals of the envelope
             and fit.
    Plot 3): Plots the V-Venv to make sure the next plot is clear as possible.
    Plot 4): Plots the integrated voltage data, fit of the integrated voltage,
             and residuals. (See DataAnalysis.py ->
             Analyzer.calculate_bunchArea_diff for the definition of
             'integrated bunch area')
    Plot 5): Plots the spectrum and shows the measured revolution frequency
             and analyzed revolution frequency (difference is that measured
             was collected during the experiment and analyzed is via the
             Fourier Transform).

    assumes:
    - analyzed_waveform is initialized properly and all of the analysis
      methods have been called, in the correct order (see DataAnalysis.py) for
      the correct order.

    - save_path, if given as an argument is a valid directory

    arguments:
    - analyzed_waveform(Analyzer): contains all of the data from the file,
                                   and analysis.
    - save_path(string): the directory where the plots are saved, defaults
                         to the current working directory.
    '''
    # Get the filename from the waveform
    filename = analyzed_waveform.waveform.filename

    ############ FIRST PLOT ############
    # Assign raw and cleaned data from the analyzed_waveform object
    t_raw, V_raw = (analyzed_waveform.waveform.raw_t,
                    analyzed_waveform.waveform.raw_V)
    t_clean, V_clean = (analyzed_waveform.waveform.clean_t,
                        analyzed_waveform.waveform.clean_V)

    # Initialize figure and axis objects
    fig = plt.figure(figsize = (10,8))
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212, sharex=ax1)

    # Plot the the raw data
    ax1.plot(t_raw, V_raw, c=color_dict['data'],
            label= 'RAW '+filename.replace('_', '\_') , alpha = 0.7)

    # Plot the cleaned data
    ax2.plot(t_clean, V_clean, c=color_dict['data'],
            label='CLEAN ' + filename.replace('_', '\_'), alpha = 0.7)

    # Begin formatting plot:
    # Change the x and y tick label frequency
    for label in ax1.get_xticklabels()[1::2]:
        label.set_visible(False)
    for label in ax1.get_yticklabels()[1::2]:
        label.set_visible(False)

    # Since ax1 and ax2 are sharing the x-axis, set visibility to false for ax1
    plt.setp(ax1.get_xticklabels(), visible=False)

    # Adding horizontal grid lines
    for ax in (ax1, ax2):
        ax.yaxis.grid(True)

        # Remove axis spines
        ax.spines["top"].set_visible(True)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(True)
        ax.spines["left"].set_visible(False)

        # Labels
        if ax != ax1:
            ax.set_xlabel('Time $[s]$', fontsize=15)
        ax.set_ylabel('Voltage $[V]$', fontsize=15)
        ax.legend(loc='best', fancybox=True, fontsize=15)

    plt.tight_layout()
    plt.savefig(save_path+'raw_and_cleaned.pdf')

    ############ SECOND PLOT ############

    envelope_t, envelope_V = (analyzed_waveform.envelope_t,
                              analyzed_waveform.envelope_V)
    envelope_t_fit, envelope_V_fit = (analyzed_waveform.envelope_t_predict,
                                      analyzed_waveform.envelope_V_predict)

    # Initialize figure and axis objects
    fig = plt.figure(figsize = (10,8))

    # Make axis 1 three times as large as axis 2
    gs = gridspec.GridSpec(2, 1, height_ratios=[3,1])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1], sharex=ax1)

   # Plot the the cleaned data
    ax1.plot(t_clean, V_clean, c=color_dict['data'],
            label='CLEAN '+filename.replace('_', '\_'), alpha = 0.7)

    # Plot the raw envelope
    ax1.plot(envelope_t, envelope_V, c='black', linestyle='None',
             marker='d', markersize=4 ,label='Envelope data', alpha = 0.7)

    # Plot the envelope fit
    ax1.plot(envelope_t_fit, envelope_V_fit, c=color_dict['fit'],
             linestyle='-', linewidth=3,
             label='Envelope fit', alpha = 0.7)

    # Plot the residuals
    ax2.plot(envelope_t_fit, np.log10(np.abs(envelope_V_fit-envelope_V)),
             c=color_dict['residuals'],
             label='$\log_{10}$($|$env-env$_{fit}$$|$)')

    # Begin formatting plot:
    # Change the x and y tick label frequency
    for label in ax1.get_xticklabels()[1::2]:
        label.set_visible(False)
    for label in ax1.get_yticklabels()[1::2]:
        label.set_visible(False)

    # Since ax1 and ax2 are sharing the x-axis, set visibility to false for ax1
    plt.setp(ax1.get_xticklabels(), visible=False)

    # Adding horizontal grid lines
    for ax in (ax1, ax2):
        ax.yaxis.grid(True)

        # Remove axis spines
        ax.spines["top"].set_visible(True)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(True)
        ax.spines["left"].set_visible(False)

    # Labels for ax1; No xlabel since it's being shared
    ax1.set_ylabel('Voltage $[V]$', fontsize=15)
    ax1.legend(loc='best', fancybox=True, fontsize=15, numpoints=1)

    # Labels for ax2
    ax2.set_xlabel('Time $[s]$')
    ax2.set_ylabel('Residuals', fontsize=15)
    ax2.legend(loc='best', fancybox=True, fontsize=15, numpoints=1)

    # Tight layout
    plt.tight_layout()

    # Save the figure
    plt.savefig(save_path+'envelope_and_residuals.pdf')

    ############ THIRD PLOT ############

    fig = plt.figure(figsize=(7,4))
    ax = plt.gca()
    
    # Plot the ideal data
    ax.plot(analyzed_waveform.ideal_t, analyzed_waveform.ideal_V,
            c = color_dict['data'],
            label='IDEAL '+filename.replace('_', '\_'), alpha = 0.7)
    
    # Begin formatting plot

    # Change xtick and ytick frequency
    for label in ax.get_xticklabels()[1::2]:
        label.set_visible(False)
    for label in ax.get_yticklabels()[1::2]:
        label.set_visible(False)

    # Adding horizontal grid lines
        ax.yaxis.grid(True)

    # Remove axis spines
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(True)
    ax.spines["left"].set_visible(False)

    # Labels
    plt.xlabel('Time $[s]$', fontsize=15)
    plt.ylabel('Voltage $[V]$', fontsize=15)
    plt.legend(loc='best', fancybox=True, fontsize=15)

    # Tight layout
    plt.tight_layout()

    # Save the plot
    plt.savefig(save_path+'ideal_v.pdf')
    
    ############ FOURTH PLOT ############

    bunch_t, bunch_V = (np.asarray(analyzed_waveform.int_bunch_diff_x)/analyzed_waveform.waveform.true_freq,
                              analyzed_waveform.int_bunch_diff)
    bunch_t_fit, bunch_V_fit = (np.asarray(analyzed_waveform.int_bunch_diff_x_predict)/analyzed_waveform.waveform.true_freq,
                                      analyzed_waveform.int_bunch_diff_predict)

    # Initialize figure and axis objects
    fig = plt.figure(figsize = (10,8))

    # Make axis 1 three times as large as axis 2
    gs = gridspec.GridSpec(2, 1, height_ratios=[3,1])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1], sharex=ax1)

    # Plot the raw int bunch diff array
    ax1.plot(bunch_t, bunch_V, c='black', linestyle='None',
             marker='d', markersize=4 ,
             label=filename.replace('_', '\_')+' Integrated Voltage Data',
             alpha = 0.7)

    # Plot the fitted int bunch diff array
    ax1.plot(bunch_t_fit, bunch_V_fit, c=color_dict['fit'],
             linestyle='-', linewidth=3,
             label='Fit', alpha = 0.7)

    # Plot the residuals; note that I exclude bunch_V's last point
    # when I subtract, that is because that point is an artefact from the analysis.
    # It is visible when I plot data against fit, for completeness.
    ax2.plot(bunch_t_fit, np.log10(np.abs(bunch_V[:-1]-bunch_V_fit)),
             c=color_dict['residuals'],
             label='$\log_{10}$($|$data-data$_{fit}$$|$)')

    # Begin formatting plot:
    # Change the x and y tick label frequency
    for label in ax1.get_xticklabels()[1::2]:
        label.set_visible(False)
    for label in ax1.get_yticklabels()[1::2]:
        label.set_visible(False)

    # Since ax1 and ax2 are sharing the x-axis, set visibility to false for ax1
    plt.setp(ax1.get_xticklabels(), visible=False)

    # Adding horizontal grid lines
    for ax in (ax1, ax2):
        ax.yaxis.grid(True)

        # Remove axis spines
        ax.spines["top"].set_visible(True)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(True)
        ax.spines["left"].set_visible(False)

    # Labels for ax1; No xlabel since it's being shared
    ax1.set_ylabel('Integrated Voltage $[V\cdot s]$', fontsize=15)
    ax1.legend(loc='best', fancybox=True, fontsize=15, numpoints=1)

    # Labels for ax2
    ax2.set_xlabel('Time $[s]$')
    ax2.set_ylabel('Residuals', fontsize=15)
    ax2.legend(loc='best', fancybox=True, fontsize=15, numpoints=1)

    # Tight layout
    plt.tight_layout()

    # Save the figure
    plt.savefig(save_path+'int_bunch_array_and_residuals.pdf')

    ############ FIFTH PLOT ############

    # Initialize figure and axis objects
    fig = plt.figure(figsize = (10,8))
    ax = plt.gca()

    # Plot the the data
    ax.fill_between(analyzed_waveform.fourier_f[analyzed_waveform.range_mask]/1e6, 0,
                        analyzed_waveform.fourier_V[analyzed_waveform.range_mask],
                        color=color_dict['data'], alpha=0.7)

    # Below is just a mess of code: the first plt.plot() call is a hack to
    # just get a legend entry for ax.fill_between() -- since the fill_between
    # plotting method does not have plt.legend() compatibility.
    # The second and third plt.plot() calls are the analyzed and measured
    # frequency, respectively.
    plt.plot(0,0,
                 color=color_dict['data'], linestyle = 'None', marker='s',
                 fillstyle = 'full', markersize = 10, markeredgecolor=color_dict['data'],
                 label=filename.replace('_', '\_')+' Spectrum')

    # Plot the analyzed frequency (corresponds to the peak in the fourier) plot
    plt.plot([analyzed_waveform.measured_revfreq/1e6,analyzed_waveform.measured_revfreq/1e6],
                  [0, max(analyzed_waveform.fourier_V[analyzed_waveform.range_mask])],
                  'r--', linewidth = 3,
                  label = 'Analyzed Revolution Frequency: %0.3f $MHz$'%(analyzed_waveform.measured_revfreq/1e6))

    # Plot the measured frequency -- the frequency when the experiment was
    # conducted
    plt.plot([analyzed_waveform.waveform.true_freq/1e6, analyzed_waveform.waveform.true_freq/1e6],
                  [0, max(analyzed_waveform.fourier_V[analyzed_waveform.range_mask])],
                  'g--', linewidth = 3,
                  label = 'Measured Revolution Frequency: {0} $MHz$'.format(analyzed_waveform.waveform.true_freq/1e6))

    # Begin formatting plot

    # Set the axis to show integration range
    ax.set_xlim([analyzed_waveform.waveform.true_freq*0.5/1e6,
                 analyzed_waveform.waveform.true_freq*1.5/1e6])

    # Change xtick and ytick frequency
    for label in ax.get_xticklabels()[2::2]:
        label.set_visible(False)
    for label in ax.get_yticklabels()[1::2]:
        label.set_visible(False)

    # Adding horizontal grid lines
        ax.yaxis.grid(True)

    # Remove axis spines
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(True)
    ax.spines["left"].set_visible(False)

    # Labels
    plt.xlabel('Frequency [MHz]', fontsize=15)
    plt.ylabel('$\mathcal{F}[V]$ $[Vs]$', fontsize=15)
    plt.legend(fancybox=True, fontsize=15, numpoints=1)

    # Tight layout
    plt.tight_layout()

    # Save the plot
    plt.savefig(save_path+'fourier_plot.pdf')
    return

def batch_analyze(data_path):
    '''
    Analyzes all the preferred files using the DataAnalysis.py module in the
    directory given by data_path.

    assumes:
    - the data files are comma separated files with time and voltage data in
      columns 3 and 4 respectively.

    arguments:
    - data_path(string): the relative or absolute path to the directory where
                         the files are stored (note in windows you need \\
                         instead of \ so that python can see it as a backslash)

    returns:
    - analyzed_batch(list): a list of Analyzer objects, ready to be plotted.
    '''
    analyzed_batch = []
    for filename in os.listdir(data_path):
        if 'inu' in filename:
            waveform = DataAnalysis.Waveform()
            waveform.load(data_path, filename)
            waveform.clean()
            analyzed_waveform = DataAnalysis.Analyzer(waveform)
            analyzed_waveform.calculate_envelope()
            analyzed_waveform.calculate_bunchArea_diff()
            analyzed_waveform.calculate_P()
            analyzed_batch.append(analyzed_waveform)
    return analyzed_batch

def batch_plot_results(analyzer_list, save_path = '.\\'):
    '''
    Plots the final results of the analysis - there are a total of X... TODO


    assumes:
    - analyzer_list contains Analyzer objects fully analyzed and ready to be
      plotted.

    arguments:
    - analyzer_list(list): a list of Analyzer objects

    - save_path(string): a string containing the directory to where the plots
                         are saved. (defaults to current working directory).
    returns:
    - nothing: only saves
    '''

    # Create the alpha list
    alpha = np.asarray([x.alpha for x in analyzer_list])

    # Create the P list, and then normalize
    P = [x.P for x in analyzer_list]
    P = np.asarray([x/max(P) for x in P])

    # Create the LINAC delay list
    LINACdelay = np.asarray([x.waveform.LINAC_delay for x in analyzer_list])

    # Create the pulse width list
    pulseWidth = np.asarray([x.waveform.pulse_width for x in analyzer_list])

    ##### PLOT 1 #####
    # Initialize figure and axis objects
    fig = plt.figure(figsize = (10,8))
    ax = plt.gca()

    # Plot the the data
    gridsize=30

    # Plot a heat map where x and y axes are the pulse width and
    # linac delay, and the color represents alpha
    plt.hexbin(LINACdelay, pulseWidth, C=alpha, gridsize=gridsize, bins=None)
    plt.axis([min(LINACdelay) * 0.99, max(LINACdelay) * 1.01,
              min(pulseWidth) * 0.95, max(pulseWidth) * 1.05])

    # Begin formatting plot

    ax.set_yticks(np.arange(np.floor(min(pulseWidth)),
                            np.ceil(max(pulseWidth))+1, 5.0))

    # Adding horizontal grid lines
    ax.yaxis.grid(True)

    # Remove axis spines
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(True)
    ax.spines["left"].set_visible(False)

    # Labels
    cb = plt.colorbar()
    cb.set_label('$\\alpha$ $[s^{-1}]$', rotation=270, fontsize=20, labelpad=30)
    plt.xlabel('LINAC Trigger Delay $[\mu s]$', fontsize=15)
    plt.ylabel('Pulse Width $[\mu s]$', fontsize=15)
    plt.legend(loc='best', fancybox=True, fontsize=15)

    # Tight layout
    plt.tight_layout()

    # Save the plot
    plt.savefig(save_path+'result_alpha.pdf')

    ##### PLOT 2 #####
    # Initialize figure and axis objects
    fig = plt.figure(figsize = (10,8))
    ax = plt.gca()

    # Plot the the data
    gridsize=30

    # Plot a heat map where x and y axes are the pulse width and
    # linac delay, and the color represents alpha
    plt.hexbin(LINACdelay, pulseWidth, C=P, gridsize=gridsize, bins=None)
    plt.axis([min(LINACdelay) * 0.99, max(LINACdelay) * 1.01,
              min(pulseWidth) * 0.95, max(pulseWidth) * 1.05])

    # Begin formatting plot

    ax.set_yticks(np.arange(np.floor(min(pulseWidth)),
                            np.ceil(max(pulseWidth))+1, 5.0))

    # Adding horizontal grid lines
    ax.yaxis.grid(True)

    # Remove axis spines
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(True)
    ax.spines["left"].set_visible(False)

    # Labels
    cb = plt.colorbar()
    cb.set_label('$\\mathbf{P}$ $[V]$', rotation=270, fontsize=20, labelpad=30)
    plt.xlabel('LINAC Trigger Delay $[\mu s]$', fontsize=15)
    plt.ylabel('Pulse Width $[\mu s]$', fontsize=15)
    plt.legend(loc='best', fancybox=True, fontsize=15)

    # Tight layout
    plt.tight_layout()

    # Save the plot
    plt.savefig(save_path+'result_P.pdf')


    ### Projection results
    
    alpha1 = alpha[np.asarray(LINACdelay) == 1570.0]
    alpha2 = alpha[np.asarray(LINACdelay) == 1590.0]

    P1 = P[np.asarray(LINACdelay) == 1570.0]
    P2 = P[np.asarray(LINACdelay) == 1590.0]
    
    x_1 = pulseWidth[np.asarray(LINACdelay) == 1570.0]
    x_2 = pulseWidth[np.asarray(LINACdelay) == 1590.0]
    
    ## PLOT ALPHA
     # Initialize figure and axis objects
    fig = plt.figure(figsize = (10,8))

    # Make the axis objects
    ax1 = plt.subplot(2,1,1)
    ax2 = plt.subplot(2,1,2)

    # Plot the projected alpha values
    ax1.plot(x_1, alpha1, c='black', linestyle='None',
             marker='o', markersize=15 ,
             label='Projected: Trigger Delay = 1570 $\mu s$',
             alpha = 0.7)
    
    # Plot the projected alpha values    
    ax2.plot(x_2, alpha2, c='black', linestyle='None',
             marker='o', markersize=15 ,
             label='Projected: Trigger Delay = 1590 $\mu s$',
             alpha = 0.7)
    
    # Adding horizontal grid lines
    for ax in (ax1, ax2):
        ax.yaxis.grid(True)

        # Remove axis spines
        ax.spines["top"].set_visible(True)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(True)
        ax.spines["left"].set_visible(False)

        # Labels for ax1; No xlabel since it's being shared
        ax.set_ylabel('$\\alpha$ $[s^{-1}]$', fontsize=15)
        ax.set_xlabel('Pulse Width $[\mu s]$')
        ax.legend(loc='best', fancybox=True, fontsize=15, numpoints=1)
    # Tight layout
    plt.tight_layout()

    # Save the figure
    plt.savefig(save_path+'alpha_projection.pdf')
    
    ## PLOT P
     # Initialize figure and axis objects
    fig = plt.figure(figsize = (10,8))

    # Make the axis objects
    ax1 = plt.subplot(2,1,1)
    ax2 = plt.subplot(2,1,2)

    # Plot the projected alpha values
    ax1.plot(x_1, P1, c='black', linestyle='None',
             marker='o', markersize=15 ,
             label='Projected: Trigger Delay = 1570 $\mu s$',
             alpha = 0.7)
    
    # Plot the projected alpha values    
    ax2.plot(x_2, P2, c='black', linestyle='None',
             marker='o', markersize=15 ,
             label='Projected: Trigger Delay = 1590 $\mu s$',
             alpha = 0.7)
    
    # Adding horizontal grid lines
    for ax in (ax1, ax2):
        ax.yaxis.grid(True)

        # Remove axis spines
        ax.spines["top"].set_visible(True)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(True)
        ax.spines["left"].set_visible(False)

        # Labels for ax1; No xlabel since it's being shared
        ax.set_ylabel('\\mathbf{P} $[V]$', fontsize=15)
        ax.set_xlabel('Pulse Width $[\mu s]$')
        ax.legend(loc='best', fancybox=True, fontsize=15, numpoints=1)
    # Tight layout
    plt.tight_layout()

    # Save the figure
    plt.savefig(save_path+'P_projection.pdf')


    return

def main():
    '''
    This is the main function -- it calls all the specified functions and
    methods and completes the analysis.
    '''
    # The following lines create the waveform object, load in the data
    # and then does the analysis. It is important that the analysis code
    # is called in this order -- with the exception of ...calculate_P()
    # this is the only analysis method which has no dependencies.
    # I should put exceptions / error messages in DataAnalysis.py to
    # outline this fact, but I'm too lazy at the moment.
    waveform = DataAnalysis.Waveform()
    waveform.load(data_path, filename_sample)
    waveform.clean(trim_start=15000)
    analyzed_waveform = DataAnalysis.Analyzer(waveform)
    analyzed_waveform.calculate_envelope()
    analyzed_waveform.calculate_bunchArea_diff()
    analyzed_waveform.calculate_P()

    # Call the plotting functions
#    plot_noise(data_path, filename_noise, save_path='.\\plots\\')
#    plot_sample_waveform(analyzed_waveform, save_path='.\\plots\\')

    # Do the analysis on ALL the data and plot the results
    analyzed_wave_list = batch_analyze(data_path)
    batch_plot_results(analyzed_wave_list)#, save_path='.\\plots\\')


if __name__=='__main__':
    main()
