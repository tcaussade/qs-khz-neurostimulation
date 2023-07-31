#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 13:42:16 2017

@author: leo
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sp
import scipy.io as sio
from .detect_peaks import detect_peaks


def raster(event_times_list, height=1, xlim=None, **kwargs):
    """
    Creates a raster plot
    Parameters
    ----------
    event_times_list : iterable
                       a list of event time iterables
    color : string
            color of vlines
    Returns
    -------
    ax : an axis containing the raster plot
    """
    ax = plt.gca()
    for ith, trial in enumerate(event_times_list):
        ymin = ith + 1 - height / 2.
        ymax = ith + 1 + height / 2.
        plt.vlines(trial, ymin, ymax, **kwargs)
    plt.ylim(1 - height / 2., len(event_times_list) + height / 2.)
    if xlim is not None:
        plt.xlim(xlim)
    return ax


def detect_spikes(samples, t, thresh=0, fs=None):
    if fs is None:
        fs = 1.0 / (t[1] - t[0])  # t in ms
    Nsamples = fs * 1  # fs in kHz, AP ~ 1 ms
    indices = detect_peaks(samples, mph=thresh, mpd=Nsamples)
    tspk = t[indices]
    return tspk


def firing_rate(spk, t, window_width, filter_type):
    '''spk: spike times
    t: time vector
    window_width: window width in ms
    filter_type: name of window to apply (string)'''

    dt = t[1] - t[0]
    N = int(window_width / dt)
    # divide by width in s to get fr in Hz
    h = sp.get_window(filter_type, N) / (window_width * 1e-3)
    spkbin = np.zeros(len(t))
    for i in range(len(spk)):
        ind = np.nonzero(t > spk[i])[0]
        if len(ind) == 0:
            print('Spike out of time interval: ' + str(spk[i]))
        else:
            ind = ind[0]
        spkbin[ind] = spkbin[ind] + 1
    fr = sp.convolve(spkbin, h, 'same')
    # crop outside of time window
    tfr = t[(t <= t[-1] - window_width / 2.) & (t >= t[0] + window_width / 2.)]
    fr = fr[(t <= t[-1] - window_width / 2.) & (t >= t[0] + window_width / 2.)]
    return fr, tfr


def plot_case(simvmid, simvend, simhs, fr, t, themodel, thefiber, theamps):

    fig = plt.figure()
    for j in range(0, len(theamps)):

        ax = fig.add_subplot(len(theamps), 2, 2 * j + 1)

        ax.plot(1e-3 * t, simvend[themodel, thefiber, j][0:len(t)], 'k')
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position('right')
        ax.tick_params(direction='out')
        ax.set_ylim(-80, 40)
        ax.set_yticks([0])
        ax.set_yticklabels([])
        if j != len(theamps) - 1:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel('Time (s)')
        if j == 0:
            ax.set_ylabel('       Amp\n       ' +
                          str(theamps[j]) + '*T', rotation=0)
            plt.title('Firing Activity', color='r')
        else:
            ax.set_ylabel('       ' + str(theamps[j]) + '*T', rotation=0)

        ax2 = plt.twinx(ax)
        ax2.plot(1e-3 * t, fr[j], 'r')
        ax2.yaxis.tick_left()
        ax2.yaxis.set_label_position('left')
        ax2.set_ylim(0, 105)
        if j != len(theamps) - 1:
            ax2.set_yticklabels([])
        else:
            ax2.set_yticks([0, 100])
            ax2.set_yticklabels([0, 100], color='r')
            ax2.set_ylabel('FR (Hz)', color='r')

        ax = fig.add_subplot(len(theamps), 2, 2 * j + 2)
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position('right')
        if len(simhs) == 0:
            ax.plot(1e-3 * t, np.ones([len(t)]), 'b')
        else:
            ax.plot(1e-3 * t, simhs[themodel, thefiber, j][0:len(t)], 'b')
        ax.set_ylim(-.1, 1.1)
        ax.set_yticks([0, 1])
        if j != len(theamps) - 1:
            ax.set_xticklabels([])
            ax.set_yticklabels([])

        else:
            ax.set_yticklabels([0, 1], color='b')
            ax.set_ylabel('$h_s$', color='b')
            ax.set_xlabel('Time (s)')
        if j == 0:
            plt.title('Slow Inactivation Gate', color='b')


def load_firing_data(matfile, simvmid, simvend, simhs,
                     frend, t, themodel, thefiber):
    matdata = sio.loadmat(matfile)
    if 'simvmid' in matdata:
        simvmid.append(matdata['simvmid'])  # name Im using in Matlab
    if 'simvend' in matdata:
        simvend.append(matdata['simvend'])  # name Im using in Matlab
        fr = []
        for j in range(0, np.shape(simvend)[-1]):
            spk = detect_spikes(np.squeeze(
                simvend[-1][themodel, thefiber, j][0:len(t)]), t, 0)
            fr.append(firing_rate(spk, t, 500, 'boxcar'))
            frend.append(fr)
    if 'simhs' in matdata:
        simhs.append(matdata['simhs'])  # name Im using in Matlab

    return simvmid, simvend, simhs, frend


def spikes_not_in(spkA, spkB, dt):

    notin = []
    common = []
    for spk in spkA:
        if not (spk in spkB or (min(abs(spkB - spk)) < dt)):
            notin.append(spk)
        else:
            common.append(spk)

    return notin, common
