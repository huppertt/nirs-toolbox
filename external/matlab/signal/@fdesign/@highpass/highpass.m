function this = highpass(varargin)
%HIGHPASS   Construct a HIGHPASS filter designer.
%   D = FDESIGN.HIGHPASS(SPECSTRING,VALUE1,VALUE2,...) constructs a
%   highpass filter designer D. Note that D is not the design itself, it
%   only contains the design specifications. In order to design the filter,
%   one needs to invoke the DESIGN method on D.
%   For example (more examples below):
%   D = fdesign.highpass('Fst,Fp,Ast,Ap',0.4,0.5,80,1);
%   H = design(D,'equiripple'); % H is a DFILT
%
%   SPECSTRING is a string that determines what design specifications will
%   be used. There are several possible specifications, a complete list is
%   given below.
%
%   Different specification types may have different design methods
%   available. Use DESIGNMETHODS to get a list of design methods
%   available for a given SPEC: designmethods(D).
%
%   VALUE1, VALUE2, etc. are scalars that provide the value of the
%   corresponding specification. In the example above, this means that Fst
%   = 0.4, Fp = 0.5, Ast = 80, and Ap = 1. Use get(D, 'description') for a
%   description of VALUE1, VALUE2, etc.
%
%   By default, all frequency specifications are assumed to be in
%   normalized frequency units. Moreover, all magnitude specifications are
%   assumed to be in dB.
%
%   D = FDESIGN.HIGHPASS(...,Fs) provides the sampling frequency of the
%   signal to be filtered. Fs must be specified as a scalar trailing the
%   other numerical values provided. For this case, Fs is assumed to be in
%   Hz as are all other frequency values provided. Note that you don't
%   change the specification string in this case. In the example above, if
%   the input signal is sampled at 8 kHz, we can obtain the same filter by
%   specifying the frequencies in Hz as:
%   D = fdesign.highpass('Fst,Fp,Ast,Ap',1600,2000,80,1,8000);
%   H = design(D,'equiripple');
%
%   D = FDESIGN.HIGHPASS(...,MAGUNITS) specifies the units for any magnitude
%   specification given. MAGUNITS can be one of the following: 'linear',
%   'dB', or 'squared'. If this argument is omitted, 'dB' is assumed. Note
%   that the magnitude specifications are always converted and stored in dB
%   regardless of how they were specified. If Fs is provided, MAGUNITS must
%   be provided after Fs in the input argument list.
%   
%   The full list of possible values for SPECSTRING (not case sensitive)
%   is:
%
%         'Fst,Fp,Ast,Ap' (minimum order; default)
%         'N,F3dB'
%         'Nb,Na,F3dB'
%         'N,F3dB,Ap' (*)
%         'N,F3dB,Ast' (*)
%         'N,F3dB,Fp' (*)
%         'N,F3dB,Ast,Ap' (*)
%         'N,Fc'
%         'N,Fc,Ast,Ap' 
%         'N,Fp,Ap'
%         'N,Fp,Ast,Ap'
%         'N,Fst,Ast'
%         'N,Fst,F3dB' (*)
%         'N,Fst,Fp'
%         'N,Fst,Ast,Ap' (*)
%         'N,Fst,Fp,Ap' (*)
%         'N,Fst,Fp,Ast' (*)
%         'Nb,Na,Fst,Fp' (*)
%
%  where 
%       Ap    - Passband Ripple (dB)
%       Ast   - Stopband Attenuation (dB)
%       F3dB  - 3dB Frequency
%       Fc    - Cutoff Frequency
%       Fp    - Passband Frequency
%       Fst   - Stopband Frequency
%       N     - Filter Order
%       Nb    - Numerator Order
%       Na    - Denominator Order
%
%   D = FDESIGN.HIGHPASS(Fstop, Fpass, Astop, Apass) uses the  default
%   SPECSTRING ('Fst,Fp,Ast,Ap') and sets the stopband-edge frequency,
%   passband-edge frequency, stopband attenuation, and passband ripple.
%
%   % Example #1 - Design an equiripple highpass filter. Specify ripples in
%   % linear units.
%   d  = fdesign.highpass('Fst,Fp,Ast,Ap',0.31,0.32,1e-3,1e-2,'linear');
%   Hd = design(d, 'equiripple');
%   fvtool(Hd)
%
%   % Example #2 - Design a Chebyshev Type I IIR filter with a passband
%   % ripple of 0.5 dB and a 3 dB cutoff frequency at 9600 Hz. (*)
%   Fs = 48000; % Sampling frequency of input signal
%   d  = fdesign.highpass('N,F3dB,Ap', 10, 9600, .5, Fs);
%   Hd = design(d, 'cheby1');
%   fvtool(Hd)
%
%   % Example #3 - Design an equiripple filter with a stopband edge of
%   % 0.65*pi rad/sample and a passband edge of 0.75*pi rad/sample. Shape the
%   % stopband to have a linear decay with a slope of 10 dB/rad/sample. (*)
%   d  = fdesign.highpass('N,Fst,Fp', 50, 0.65, 0.75);
%   Hd = design(d, 'equiripple','StopbandShape','linear','StopbandDecay',20);
%   fvtool(Hd)
%
%   %(*) DSP System Toolbox required
%
%   See also FDESIGN, FDESIGN/SETSPECS, FDESIGN/DESIGN, FDESIGN/DESIGNOPTS.

%   Copyright 1999-2013 The MathWorks, Inc.

this = fdesign.highpass;

[varargin,flag] = finddesignfiltflag(this,varargin);

set(this, 'Response', 'Highpass');

if flag 
  specObj = this.getcurrentspecs;
  specObj.FromDesignfilt = true;
end

this.setspecs(varargin{:});

capture(this);

% [EOF]
