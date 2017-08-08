function pwr = avgpower(this,freqrange)
%AVGPOWER    Average power.
%
%   AVGPOWER is not recommended.  Use <a href="matlab:help bandpower">bandpower</a> instead.
%
%   AVGPOWER(H) computes the average power in a signal via a rectangle
%   approximation of the integral of the Power Spectral Density (PSD) of
%   such signal, given in H. If the PSD data in H is a matrix the operation
%   will be done on each column.  H is a DSP data (<a href="matlab:help dspdata">dspdata</a>) object.
%
%   Depending on the value of the SpectrumType property of H, AVGPOWER
%   calculates the total average power either over the one-sided or
%   two-sided spectrum.  For a one-sided spectrum, the computation is done
%   over the range [0,pi] for an even number of frequency points, and over
%   the range [0,pi) for an odd number of frequency points.  For a
%   two-sided spectrum, the calculation is done over the range [0,2pi).
%
%   AVGPOWER(H,FREQRANGE) specifies the frequency range on which to measure
%   the average power.  FREQRANGE is a two-element vector of real values,
%   specifying the two frequencies between which you want to measure the
%   power.  
%
%   NOTE: If the frequency range values don't exactly match the frequency
%   values stored in the object H then the next closest value is used.
%   Also, the first element in the two element frequency range vector is
%   included in calculating the average power while the second value is
%   excluded.
%
%   EXAMPLE 
%      Calculate the average power of a sinusoid with a frequency component
%      at .2 Hz.  The average power (in one period) is given by A^2/2 = 8,
%      where A is the amplitude of the sinusoid.
%
%      f = .2; Fs = 1;
%      n = 0:999;
%      x = 4*sin(2*pi*f/Fs*n);
%
%      h = spectrum.periodogram('rectangular');     % Periodogram spectral estimator. 
%      hopts = psdopts(h);                          % Default options.
%      set(hopts,'NFFT',2^12,'Fs',Fs,'SpectrumType','twosided');
%      hpsd = psd(h,x,hopts);                       % PSD data object.
%      pwr = avgpower(hpsd)
%
%      % Calculating and plotting the one-sided PSD we notice that most of
%      % the power lies in the range of 0.1 and 0.3.  Therefore, 
%      % calculating the power over that range should return most, but not
%      % all, of the power of the signal.  This is due to spectral leakage.
%
%      hopts.SpectrumType='onesided';
%      hpsd = psd(h,x,hopts);
%      plot(hpsd);
%      pwr = avgpower(hpsd,[.1 .3])

%   Author(s): P. Pacheco
%   Copyright 1988-2012 The MathWorks, Inc.

% Help for the AVGPOWER method.

% [EOF]
