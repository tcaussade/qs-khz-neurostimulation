#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 13:42:16 2017

@author: leo
"""
import numpy as np
import matplotlib.pyplot as plt
from .detect_peaks import detect_peaks

def detect_spikes(samples, t, thresh=0, fs=None):
    if fs is None:
        fs = 1.0 / (t[1] - t[0])  # t in ms
    Nsamples = fs * 1  # fs in kHz, AP ~ 1 ms
    indices = detect_peaks(samples, mph=thresh, mpd=Nsamples)
    tspk = t[indices]
    return tspk