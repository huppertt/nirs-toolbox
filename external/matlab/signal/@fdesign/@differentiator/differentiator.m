function this = differentiator(varargin)
%DIFFERENTIATOR   Construct a differentiator filter designer.
%   D = FDESIGN.DIFFERENTIATOR constructs a differentiator filter designer D.
%
%   D = FDESIGN.DIFFERENTIATOR(SPEC) initializes the filter designer
%   'Specification' property to SPEC.  SPEC is one of the following
%   strings and is not case sensitive:
%
%       'N'             - Full band differentiator (default) 
%       'N,Fp,Fst'      - Partial band differentiator
%       'N,Fp,Fst,Ap'   - Partial band differentiator (*)
%       'N,Fp,Fst,Ast'  - Partial band differentiator (*)
%       'Ap'            - Minimum order full band differentiator (*)
%       'Fp,Fst,Ap,Ast' - Minimum order partial band differentiator (*)
%
%  where 
%       Ap    - Passband Ripple (dB)
%       Ast   - Stopband Attenuation (dB)
%       Fp    - Passband Frequency
%       Fst   - Stopband Frequency
%       N     - Filter Order
%
%   By default, all frequency specifications are assumed to be in
%   normalized frequency units. Moreover, all magnitude specifications are
%   assumed to be in dB.
% 
%   Different specification may have different design methods available.
%   Use DESIGNMETHODS(D) to get a list of design methods available for a
%   given SPEC.
%
%   D = FDESIGN.DIFFERENTIATOR(SPEC, SPEC1, SPEC2, ...) initializes the
%   filter designer specifications with SPEC1, SPEC2, etc.
%   Use GET(D,'DESCRIPTION') for a description of SPEC1, SPEC2, etc.
%
%   D = FDESIGN.DIFFERENTIATOR(N) uses the default SPEC ('N') and sets 
%   filter order. 
%
%   D = FDESIGN.DIFFERENTIATOR(...,Fs) specifies the sampling frequency
%   (in Hz). In this case, all other frequency specifications are also in
%   Hz.
%
%   D = FDESIGN.DIFFERENTIATOR(...,MAGUNITS) specifies the units for any
%   magnitude specification given in the constructor. MAGUNITS can be one
%   of the following: 'linear', 'dB', or 'squared'. If this argument is
%   omitted, 'dB' is assumed. Note that the magnitude specifications are
%   always converted and stored in dB regardless of how they were
%   specified.
%
%   % Example #1 - Design a 33rd order full band differentiator.
%   d = fdesign.differentiator(33);
%   designmethods(d);
%   Hd = design(d,'firls');
%   fvtool(Hd,'MagnitudeDisplay','Zero-phase','FrequencyRange','[-pi, pi)')
%
%   % Example #2 - Design a narrow band differentiator.
%   %              Differentiate the 25% lowest frequencies of the Nyquist 
%   %              range and filter the higher frequencies. 
%   d = fdesign.differentiator('N,Fp,Fst',54,.25,.3);
%   Hd = design(d,'equiripple'); 
%   % Weight the stopband to increase the stopband attenuation
%   Hd1 = design(d,'equiripple','Wstop',4); 
%   fvtool(Hd,Hd1,'MagnitudeDisplay','Zero-phase','FrequencyRange','[-pi, pi)','legend','on')
%
%   % Example #3 - Design a minimum order wide band differentiator. (*)
%   d = fdesign.differentiator('Fp,Fst,Ap,Ast',.8,.9,1,80);
%   designmethods(d);
%   Hd = design(d,'equiripple'); 
%   fvtool(Hd,'MagnitudeDisplay','Zero-phase','FrequencyRange','[-pi, pi)')
%
%   %(*) DSP System Toolbox required
%
%   See also FDESIGN, FDESIGN/SETSPECS, FDESIGN/DESIGN.

%   Copyright 2005-2013 The MathWorks, Inc.

this = fdesign.differentiator;

[varargin,flag] = finddesignfiltflag(this,varargin);

if flag 
  specObj = this.getcurrentspecs;
  specObj.FromDesignfilt = true;
end

set(this, 'Response', 'Differentiator');

this.setspecs(varargin{:});

capture(this);

% [EOF]
