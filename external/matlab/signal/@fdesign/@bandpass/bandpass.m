function this = bandpass(varargin)
%BANDPASS   Construct a BANDPASS filter designer.
%   D = FDESIGN.BANDPASS(SPECSTRING,VALUE1,VALUE2,...) constructs a
%   bandpass filter designer D. Note that D is not the design itself, it
%   only contains the design specifications. In order to design the filter,
%   one needs to invoke the DESIGN method on D.
%   For example (more examples below):
%   D = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',.4,.5,.6,.7,60,1,80);
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
%   corresponding specification. Use get(D, 'description') for a
%   description of VALUE1, VALUE2, etc.
%
%   By default, all frequency specifications are assumed to be in
%   normalized frequency units. Moreover, all magnitude specifications are
%   assumed to be in dB.
%
%   D = FDESIGN.BANDPASS(...,Fs) provides the sampling frequency of the
%   signal to be filtered. Fs must be specified as a scalar trailing the
%   other numerical values provided. For this case, Fs is assumed to be in
%   Hz as are all other frequency values provided. Note that you don't
%   change the specification string in this case. In the example above, if
%   the input signal is sampled at 8 kHz, we can obtain the same filter by
%   specifying the frequencies in Hz as:
%   D = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',...
%          1600,2000,2400,2800,60,1,80,8000);
%   H = design(D,'equiripple'); 
%
%   D = FDESIGN.BANDPASS(...,MAGUNITS) specifies the units for any
%   magnitude specification given. MAGUNITS can be one of the following:
%   'linear', 'dB', or 'squared'. If this argument is omitted, 'dB' is
%   assumed. Note that the magnitude specifications are always converted
%   and stored in dB regardless of how they were specified. If Fs is
%   provided, MAGUNITS must be provided after Fs in the input argument
%   list.
%   
%   The full list of possible values for SPECSTRING (not case sensitive)
%   is:
%      'Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2' (minimum order; default)
%      'N,F3dB1,F3dB2'       
%      'N,F3dB1,F3dB2,Ap' (*)              
%      'N,F3dB1,F3dB2,Ast' (*)
%      'N,F3dB1,F3dB2,BWp' (*)
%      'N,F3dB1,F3dB2,BWst' (*)
%      'N,F3dB1,F3dB2,Ast1,Ap,Ast2' (*)
%      'N,Fc1,Fc2'
%      'N,Fc1,Fc2,Ast1,Ap,Ast2' 
%      'N,Fp1,Fp2,Ap'
%      'N,Fp1,Fp2,Ast1,Ap,Ast2' 
%      'N,Fst1,Fp1,Fp2,Fst2'
%      'N,Fst1,Fp1,Fp2,Fst2,C' (*)
%      'N,Fst1,Fst2,Ast'
%      'N,Fst1,Fp1,Fp2,Fst2,Ap' (*)
%      'Nb,Na,Fst1,Fp1,Fp2,Fst2' (*)
%
%  where 
%      Ap    - Passband Ripple (dB)
%      Ast   - Stopbands Attenuation (dB)
%      Ast1  - First Stopband Attenuation (dB)
%      Ast2  - Second Stopband Attenuation (dB)
%      BWp   - Passband Frequency Width
%      BWst  - Stopband Frequency Width
%      F3dB1 - First 3dB Frequency
%      F3dB2 - Second 3dB Frequency
%      Fc1   - First Cutoff Frequency
%      Fc2   - Second Cutoff Frequency
%      Fp1   - First Passband Frequency
%      Fp2   - Second Passband Frequency
%      Fst1  - First Stopband Frequency
%      Fst2  - Second Stopband Frequency
%      N     - Filter Order
%      Nb    - Numerator Order
%      Na    - Denominator Order
%      C     - Constrained Band Flag
%
%   D = FDESIGN.BANDPASS(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass,
%   Astop2) uses the default SPEC ('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2')
%   and sets the lower stopband-edge frequency, the lower passband-edge
%   frequency, the upper passband-edge frequency, the upper stopband-edge
%   frequency, the lower stopband attenuation, the passband ripple, and the
%   upper stopband attenuation.
%
%   % Example #1 - Design a minimum order elliptic bandpass filter.
%   d = fdesign.bandpass(.3, .4, .6, .7, 80, .5, 60);
%   Hd = design(d, 'ellip');
%   info(Hd)
%
%   % Example #2 - Design an FIR least-squares bandpass filter.
%   d = fdesign.bandpass('N,Fst1,Fp1,Fp2,Fst2',100,.3, .4, .6, .7);
%   Hd = design(d, 'firls','wstop1',100,'wpass',1,'wstop2',1);
%   fvtool(Hd)
%
%   % Example #3 - Specify frequencies in Hz.
%   d = fdesign.bandpass('N,Fp1,Fp2,Ap', 10, 9600, 14400, .5, 48000);
%   designmethods(d);
%   Hd = design(d, 'cheby1');
%
%   % Example #4 - Specify the magnitude specifications in squared units
%   d = fdesign.bandpass(.4, .5, .6, .7, .02, .98, .01, 'squared');
%   Hd = design(d, 'cheby2');
%   fvtool(Hd,'MagnitudeDisplay','Magnitude Squared');
%
%   % Example #5 - Design a constrained band equiripple bandpass filter (*)
%   d = fdesign.bandpass('N,Fst1,Fp1,Fp2,Fst2,C',60,0.3,0.4,0.7,0.8);
%   % Set constraints for the attenuation in the stopbands
%   d.Stopband1Constrained = true;
%   d.Astop1 = 60;
%   d.Stopband2Constrained = true;
%   d.Astop2 = 70;
%   Hd = design(d,'equiripple');
%   fvtool(Hd);
%
%  % (*) DSP System Toolbox required
%
%   See also FDESIGN, FDESIGN/SETSPECS, FDESIGN/DESIGN, FDESIGN/DESIGNOPTS.

%   Copyright 1999-2013 The MathWorks, Inc.

this = fdesign.bandpass;

[varargin,flag] = finddesignfiltflag(this,varargin);

set(this, 'Response', 'Bandpass');

if flag 
  specObj = this.getcurrentspecs;
  specObj.FromDesignfilt = true;
end

this.setspecs(varargin{:});

capture(this);

% [EOF]
