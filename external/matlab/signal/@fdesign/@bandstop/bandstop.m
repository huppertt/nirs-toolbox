function this = bandstop(varargin)
%BANDSTOP Construct a BANDSTOP filter designer.
%   D = FDESIGN.BANDSTOP(SPECSTRING,VALUE1,VALUE2,...) constructs a
%   bandstop filter designer D. Note that D is not the design itself, it
%   only contains the design specifications. In order to design the filter,
%   one needs to invoke the DESIGN method on D.
%   For example (more examples below):
%   D = fdesign.bandstop('Fp1,Fst1,Fst2,Fp2,Ap1,Ast,Ap2',.4,.5,.6,.7,1,80,.5);
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
%   D = FDESIGN.BANDSTOP(...,Fs) provides the sampling frequency of the
%   signal to be filtered. Fs must be specified as a scalar trailing the
%   other numerical values provided. For this case, Fs is assumed to be in
%   Hz as are all other frequency values provided. Note that you don't
%   change the specification string in this case. In the example above, if
%   the input signal is sampled at 8 kHz, we can obtain the same filter by
%   specifying the frequencies in Hz as:
%   D = fdesign.bandstop('Fp1,Fst1,Fst2,Fp2,Ap1,Ast,Ap2',...
%          1600,2000,2400,2800,1,80,.5,8000);
%   H = design(D,'equiripple'); 
%
%   D = FDESIGN.BANDSTOP(...,MAGUNITS) specifies the units for any
%   magnitude specification given. MAGUNITS can be one of the following:
%   'linear', 'dB', or 'squared'. If this argument is omitted, 'dB' is
%   assumed. Note that the magnitude specifications are always converted
%   and stored in dB regardless of how they were specified. If Fs is
%   provided, MAGUNITS must be provided after Fs in the input argument
%   list.
%   
%   The full list of possible values for SPECSTRING (not case sensitive)
%   is:
%      'Fp1,Fst1,Fst2,Fp2,Ap1,Ast,Ap2' (minimum order; default)
%      'N,F3dB1,F3dB2'
%      'N,F3dB1,F3dB2,Ap' (*)
%      'N,F3dB1,F3dB2,Ast' (*)
%      'N,F3dB1,F3dB2,BWp' (*)
%      'N,F3dB1,F3dB2,BWst' (*)
%      'N,F3dB1,F3dB2,Ap,Ast' (*)
%      'N,Fc1,Fc2'
%      'N,Fc1,Fc2,Ap1,Ast,Ap2'
%      'N,Fp1,Fp2,Ap'
%      'N,Fp1,Fp2,Ap,Ast'
%      'N,Fp1,Fst1,Fst2,Fp2'
%      'N,Fp1,Fst1,Fst2,Fp2,C' (*)
%      'N,Fp1,Fst1,Fst2,Fp2,Ap' (*)
%      'N,Fst1,Fst2,Ast'
%      'Nb,Na,Fp1,Fst1,Fst2,Fp2' (*)
%
%  where 
%      Ap    - Passbands Ripple (dB)
%      Ap1   - First Passband Ripple (dB)
%      Ap2   - Second Passband Ripple (dB)
%      Ast   - Stopband Attenuation (dB)
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
%   D = FDESIGN.BANDSTOP(Fpass1, Fstop1, Fstop2, Fpass2, Apass1, Astop,
%   Apass2) uses the default SPEC ('Fp1,Fst1,Fst2,Fp2,Ap1,Ast,Ap2')
%   and sets the lower stopband-edge frequency, the lower passband-edge
%   frequency, the upper passband-edge frequency, the upper stopband-edge
%   frequency, the lower stopband attenuation, the passband ripple, and the
%   upper passband attenuation.
%
%   % Example #1 - Design a minimum order elliptic bandstop filter.
%   d = fdesign.bandstop(.3, .4, .6, .7, .5, 60, 1);
%   designmethods(d);
%   Hd = design(d, 'ellip');
%   fvtool(Hd)
%
%   % Example #2 - Design an FIR least-squares bandstop filter.
%   d = fdesign.bandstop('N,Fp1,Fst1,Fst2,Fp2',40,.3, .4, .6, .7);
%   Hd = design(d, 'firls','wpass1',1000,'wstop',1,'wpass2',1);
%   fvtool(Hd)
%
%   % Example #3 - Specify frequencies in Hz.
%   d = fdesign.bandstop('N,Fp1,Fp2,Ap', 10, 9600, 14400, .5, 48000);
%   designmethods(d);
%   design(d, 'cheby1');
%
%   % Example #4 - Specify the magnitude specifications in squared units
%   d = fdesign.bandstop(.4, .5, .6, .7, .98, .01, .99, 'squared');
%   Hd = design(d, 'cheby2');
%   fvtool(Hd,'MagnitudeDisplay','Magnitude Squared');
%
%   % Example #5 - Design a constrained band equiripple bandstop filter (*)
%   d = fdesign.bandstop('N,Fp1,Fst1,Fst2,Fp2,C',60,0.3,0.4,0.7,0.8);
%   % Set constraints for the ripple in the passbands
%   d.Passband1Constrained = true;
%   d.Apass1 = 0.1;
%   d.Passband2Constrained = true;
%   d.Apass2 = 0.5;
%   Hd = design(d,'equiripple');
%   fvtool(Hd);
%
%   %(*) DSP System Toolbox required
%
%   See also FDESIGN, FDESIGN/SETSPECS, FDESIGN/DESIGN, FDESIGN/DESIGNOPTS.

%   Copyright 1999-2013 The MathWorks, Inc.

this = fdesign.bandstop;

[varargin,flag] = finddesignfiltflag(this,varargin);

set(this, 'Response', 'Bandstop');

if flag 
  specObj = this.getcurrentspecs;
  specObj.FromDesignfilt = true;
end

this.setspecs(varargin{:});

capture(this);

% [EOF]
