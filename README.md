# Polyphase Filter Banks

## Author

Rishav Karki  

## Overview

This project explores multiple methods for implementing 8-to-1 channelization filters for a 7-channel QAM signal using FIR and IIR filtering techniques. Each method aims to digitally downconvert each channel to baseband, apply low-pass filtering, and perform 8-to-1 downsampling to obtain outputs at 2 samples per symbol.

## Implemented Options

### Option 1: FIR Filtering with Heterodyne Downconversion

- Used 7 complex heterodynes to downconvert each QAM channel to baseband.
- Applied a 64-tap low-pass FIR filter (designed using the Remez algorithm) to each baseband signal.
- Downsampled each filtered signal by a factor of 8.

### Option 2: 8-Path Polyphase FIR Filter

- Designed a prototype 64-tap low-pass FIR filter using the Remez algorithm.
- Partitioned the filter into 8 polyphase components to perform 8-to-1 downsampling.
- Used an 8-point IFFT on the filter outputs to separate the 7 aliased channels.

### Option 3: IIR Filtering with Heterodyne Downconversion

- Used 7 complex heterodynes to downconvert each QAM channel to baseband.
- Applied a linear phase 8-path recursive IIR filter (each path has 3-coefficient all-pass filters with 24 delays).
- Did not implement Noble identity in this option.

### Option 4: 8-Path Polyphase Linear Phase IIR Filter

- Designed an 8-path IIR filter using the `lineardesign_2` script with 3 all-pass coefficients and 24 delays per path.
- Implemented 8-to-1 downsampling using the Noble identity.
- Applied an IFFT on the 8 paths to obtain outputs for the aliased channels.

## Computational Load Comparison

- **Option 1:** 65 multiplications per input sample per channel (least efficient).
- **Option 2:** 1.57 multiplications per input sample per channel.
- **Option 3:** 22 multiplications per input sample per channel.
- **Option 4:** 0.81 multiplications per input sample per channel (most efficient).

## Summary

All four filtering approaches were implemented and tested on a 7-channel QAM input signal. Frequency and time-domain responses were evaluated using MATLAB plots. Option 4, the 8-path polyphase linear phase IIR filter, proved to be the most computationally efficient, while Option 1 had the highest processing load. Each method was verified with impulse responses, downsampled outputs, and channel isolation performance.

## Files

- MATLAB scripts implementing all four options.
- Filter design coefficients and frequency response plots.
- Simulation results for impulse and QAM signal inputs.

## Tools Used

- MATLAB

