function this = isinclp(varargin)
%ISINCLP   Construct an inverse-sinc lowpass filter designer.
%   D = FDESIGN.ISINCLP constructs an inverse-sinc lowpass filter designer D.
%
%   D = FDESIGN.ISINCLP(SPEC) initializes the filter designer
%   'Specification' property to SPEC.  SPEC is one of the following
%   strings and is not case sensitive:
%
%       'Fp,Fst,Ap,Ast' - (minimum order, default) 
%       'N,Fc,Ap,Ast'   
%       'N,Fp,Ap,Ast'   
%       'N,Fp,Fst'  
%       'N,Fst,Ap,Ast'  
%
%  where 
%       Ap    - Passband Ripple (dB)
%       Ast   - Stopband Attenuation (dB)
%       Fc    - Cutoff Frequency
%       Fp    - Passband Frequency
%       Fst   - Stopband Frequency
%       N     - Filter Order
%
%   By default, all frequency specifications are assumed to be in
%   normalized frequency units. Moreover, all magnitude specifications are
%   assumed to be in dB.
%
%   Different specification types may have different design methods
%   available. Use DESIGNMETHODS(D) to get a list of design methods
%   available for a given SPEC.
%
%   D = FDESIGN.ISINCLP(SPEC, SPEC1, SPEC2, ...) initializes the filter
%   designer specifications with SPEC1, SPEC2, etc...
%   Use GET(D, 'DESCRIPTION') for a description of SPEC1, SPEC2, etc.
%
%   D = FDESIGN.ISINCLP(Fpass, Fstop, Apass, Astop) uses the default
%   SPEC ('Fp,Fst,Ap,Ast') and sets the passband-edge frequency,
%   stopband-edge frequency, passband ripple, and stopband attenuation.
%
%   D = FDESIGN.ISINCLP(..., Fs) specifies the sampling frequency (in Hz).
%   In this case, all other frequency specifications are also in Hz.
%
%   D = FDESIGN.ISINCLP(..., MAGUNITS) specifies the units for any
%   magnitude specification given in the constructor. MAGUNITS can be one
%   of the following: 'linear', 'dB', or 'squared'. If this argument is
%   omitted, 'dB' is assumed.  Note that the magnitude specifications are
%   always converted and stored in dB regardless of how they were
%   specified.
%
%   The design method of FDESIGN.ISINCLP implements a filter with a
%   passband magnitude response equal to H(w) = sinc(C*w)^(-P). You can
%   control the values of the sinc frequency factor, C, and the sinc power,
%   P, using the 'SincFrequencyFactor' and 'SincPower' options in the
%   design method. C and P default to 0.5 and 1 respectively.
%
%   % Example #1 - Design a minimum order inverse-sinc lowpass filter.
%   d = fdesign.isinclp('Fp,Fst,Ap,Ast',.4,.5,.01,40);
%   Hd = design(d);
%   fvtool(Hd, 'MagnitudeDisplay', 'Magnitude');
%
%   % Example #2 - Design a 50th order inverse-sinc lowpass filter. Set the
%   %              sinc frequency factor to 0.25, and the sinc power to 16
%   %              to achieve a magnitude response in the passband of the form 
%   %              H(w) = sinc(0.25*w)^(-16).
%   d = fdesign.isinclp('N,Fp,Fst',50,.4,.5);
%   Hd = design(d,'SincFrequencyFactor',0.25,'SincPower',16);
%   fvtool(Hd, 'MagnitudeDisplay', 'Magnitude');
%
%
%   See also FDESIGN, FDESIGN/SETSPECS, FDESIGN/DESIGN.

%   Copyright 2005-2011 The MathWorks, Inc.

this = fdesign.isinclp;

set(this, 'Response', 'Inverse-sinc Lowpass');

setspecs(this, varargin{:});

% [EOF]
