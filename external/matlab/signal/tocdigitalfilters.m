function tocdigitalfilters
% Discrete-time filter design, analysis, and implementation           
% ---------------------------------------------------------
% Specification-based filter design
%   <a href="matlab:help designfilt">designfilt</a> - Design a digital filter based on a set of specifications
%
% FIR filter design functions
%   <a href="matlab:help cfirpm">cfirpm</a>     - Complex and nonlinear phase equiripple FIR filter design
%   <a href="matlab:help fir1">fir1</a>       - Window based FIR filter design - low, high, band, stop, multi
%   <a href="matlab:help fir2">fir2</a>       - FIR arbitrary shape filter design using the frequency sampling method
%   <a href="matlab:help fircls">fircls</a>     - Constrained Least Squares filter design - arbitrary response
%   <a href="matlab:help fircls1">fircls1</a>    - Constrained Least Squares FIR filter design - low and highpass
%   <a href="matlab:help firls">firls</a>      - Optimal least-squares FIR filter design
%   <a href="matlab:help firpm">firpm</a>      - Parks-McClellan optimal equiripple FIR filter design
%   <a href="matlab:help firpmord">firpmord</a>   - Parks-McClellan optimal equiripple FIR order estimator
%   <a href="matlab:help intfilt">intfilt</a>    - Interpolation FIR filter design
%   <a href="matlab:help kaiserord">kaiserord</a>  - Kaiser window design based filter order estimation
%   <a href="matlab:help sgolay">sgolay</a>     - Savitzky-Golay FIR smoothing filter design
%
% Communications filter design functions
%   <a href="matlab:help rcosdesign">rcosdesign</a>    - Raised cosine FIR filter design
%   <a href="matlab:help gaussdesign">gaussdesign</a>   - Gaussian FIR Pulse-Shaping Filter Design
%
% IIR filter design functions
%   <a href="matlab:help butter">butter</a>     - Butterworth filter design
%   <a href="matlab:help cheby1">cheby1</a>     - Chebyshev Type I filter design (passband ripple)
%   <a href="matlab:help cheby2">cheby2</a>     - Chebyshev Type II filter design (stopband ripple)
%   <a href="matlab:help ellip">ellip</a>      - Elliptic filter design
%   <a href="matlab:help maxflat">maxflat</a>    - Generalized Butterworth lowpass filter design
%   <a href="matlab:help yulewalk">yulewalk</a>   - Yule-Walker filter design
%
% IIR filter order estimation
%   <a href="matlab:help buttord">buttord</a>    - Butterworth filter order estimation
%   <a href="matlab:help cheb1ord">cheb1ord</a>   - Chebyshev Type I filter order estimation
%   <a href="matlab:help cheb2ord">cheb2ord</a>   - Chebyshev Type II filter order estimation
%   <a href="matlab:help ellipord">ellipord</a>   - Elliptic filter order estimation
%
% Filter analysis
%   <a href="matlab:help abs">abs</a>        - Magnitude
%   <a href="matlab:help angle">angle</a>      - Phase angle
%   <a href="matlab:help filt2block">filt2block</a> - Generate Simulink filter block
%   <a href="matlab:help filternorm">filternorm</a> - 2-norm or inf-norm of a digital filter
%   <a href="matlab:help filtord">filtord</a>    - Filter order
%   <a href="matlab:help firtype">firtype</a>    - Type of linear phase FIR filter
%   <a href="matlab:help freqz">freqz</a>      - Z-transform frequency response
%   <a href="matlab:help fvtool">fvtool</a>     - Filter Visualization Tool
%   <a href="matlab:help grpdelay">grpdelay</a>   - Group delay 
%   <a href="matlab:help impz">impz</a>       - Impulse response 
%   <a href="matlab:help isallpass">isallpass</a>  - True for all-pass filters
%   <a href="matlab:help islinphase">islinphase</a> - True for linear phase filters
%   <a href="matlab:help ismaxphase">ismaxphase</a> - True for maximum phase filters
%   <a href="matlab:help isminphase">isminphase</a> - True for minimum phase filters
%   <a href="matlab:help isstable">isstable</a>   - True for stable filters
%   <a href="matlab:help phasedelay">phasedelay</a> - Phase delay 
%   <a href="matlab:help phasez">phasez</a>     - Phase response 
%   <a href="matlab:help stepz">stepz</a>      - Step response 
%   <a href="matlab:help unwrap">unwrap</a>     - Unwrap phase angle
%   <a href="matlab:help zerophase">zerophase</a>  - Zero-phase response of real filter
%   <a href="matlab:help zplane">zplane</a>     - Discrete pole-zero plot
%
% Filter implementation
%   <a href="matlab:help conv">conv</a>       - Convolution
%   <a href="matlab:help conv2">conv2</a>      - 2-D convolution
%   <a href="matlab:help convmtx">convmtx</a>    - Convolution matrix
%   <a href="matlab:help deconv">deconv</a>     - Deconvolution
%   <a href="matlab:help fftfilt">fftfilt</a>    - Overlap-add filter implementation
%   <a href="matlab:help filter">filter</a>     - Filter implementation
%   <a href="matlab:help filter2">filter2</a>    - Two-dimensional digital filtering
%   <a href="matlab:help filtfilt">filtfilt</a>   - Zero-phase version of filter
%   <a href="matlab:help filtic">filtic</a>     - Determine filter initial conditions
%   <a href="matlab:help latcfilt">latcfilt</a>   - Lattice filter implementation
%   <a href="matlab:help medfilt1">medfilt1</a>   - 1-Dimensional median filtering
%   <a href="matlab:help sgolayfilt">sgolayfilt</a> - Savitzky-Golay filter implementation
%   <a href="matlab:help sosfilt">sosfilt</a>    - Second-order sections (biquad) filter implementation
%   <a href="matlab:help upfirdn">upfirdn</a>    - Upsample, FIR filter, downsample
%
% <a href="matlab:help signal">Signal Processing Toolbox TOC</a>

%   Copyright 2005-2013 The MathWorks, Inc.


