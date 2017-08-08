function noisepsd(this)
%NOISEPSD   Power spectral density of filter output due to roundoff noise.
%   Hpsd = NOISEPSD(Hd,L) computes the power spectral density (PSD) at the
%   output of filter Hd due to roundoff noise produced by quantization
%   errors within the filter. L is the number of trials used. The PSD is
%   computed from an average over the L trials. The more trials the better
%   the estimate (at the expense of longer computation). If L is not
%   specified, it defaults to 10 trials.
%
%   Hpsd is a PSD data object. The PSD vector can be extracted from Hpsd by
%   GET(Hpsd,'Data'). It can be plotted using PLOT(Hpsd). The average power
%   of the output noise (the integral of the PSD) can be computed with
%   AVGPOWER(Hpsd).
%
%   Hpsd = NOISEPSD(Hd,L,P1,V1,P2,V2,...) specifies optional parameters via
%   parameter-value pairs. Valid pairs are:
%        Parameter           Default        Description/Valid values
%   ---------------------  -----------  ----------------------------------
%   'NFFT'                 512          Number of FFT points.
%   'NormalizedFrequency'  true         {true,false}
%   'Fs'                   'Normalized' Sampling frequency. Only used when
%                                       'NormalizedFrequency' is false.
%   'SpectrumType'         'Onesided'   {'Onesided','Twosided'}
%   'CenterDC'             false        {true,false} 
%
%   NOISEPSD(Hd,L,OPTS) uses an options object to specify the optional
%   parameters in lieu of specifying parameter-value pairs. The OPTS object
%   can be created with OPTS = NOISEPSDOPTS(Hd). Settings can be changed in
%   OPTS before calling NOISEPSD, e.g. set(OPTS,'Fs',48e3);
%
%   % Example: Compute the PSD of the output noise due to quantization in a
%   % fixed-point FIR filter implemented in direct-form.
%   f = fdesign.lowpass('N,Fp,Fst',27,.4,.6);
%   h = design(f,'equiripple','FilterStructure','dffir');
%   h.Arithmetic = 'fixed';
%   specifyall(h); h.OutputFracLength=15;
%   Hpsd = noisepsd(h);
%   plot(Hpsd)
%
%   See also DFILT/NOISEPSDOPTS, DFILT/FREQRESPEST, DFILT/SCALE,
%   DFILT/FUNCTIONS, DSPDATA/PSD.

%   Author(s): R. Losada
%   Copyright 1988-2006 The MathWorks, Inc.



% [EOF]
