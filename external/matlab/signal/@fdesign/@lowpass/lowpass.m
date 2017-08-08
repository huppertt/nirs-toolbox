function this = lowpass(varargin)
%LOWPASS   Construct a lowpass filter designer.
%   D = FDESIGN.LOWPASS(SPECSTRING,VALUE1,VALUE2,...) constructs a lowpass
%   filter designer D. Note that D is not the design itself, it only
%   contains the design specifications. In order to design the filter, one
%   needs to invoke the DESIGN method on D.
%   For example (more examples below):
%   D = fdesign.lowpass('Fp,Fst,Ap,Ast',0.4,0.5,1,80);
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
%   corresponding specification. In the example above, this means that Fp =
%   0.4, Fst = 0.5, Ap = 1, and Ast = 80. Use get(D, 'description') for a
%   description of VALUE1, VALUE2, etc.
%
%   By default, all frequency specifications are assumed to be in
%   normalized frequency units. Moreover, all magnitude specifications are
%   assumed to be in dB.
%
%   D = FDESIGN.LOWPASS(...,Fs) provides the sampling frequency of the
%   signal to be filtered. Fs must be specified as a scalar trailing the
%   other numerical values provided. For this case, Fs is assumed to be in
%   Hz as are all other frequency values provided. Note that you don't
%   change the specification string in this case. In the example above, if
%   the input signal is sampled at 8 kHz, we can obtain the same filter by
%   specifying the frequencies in Hz as:
%   D = fdesign.lowpass('Fp,Fst,Ap,Ast',1600,2000,1,80,8000);
%   H = design(D,'equiripple');
%
%   D = FDESIGN.LOWPASS(...,MAGUNITS) specifies the units for any magnitude
%   specification given. MAGUNITS can be one of the following: 'linear',
%   'dB', or 'squared'. If this argument is omitted, 'dB' is assumed. Note
%   that the magnitude specifications are always converted and stored in dB
%   regardless of how they were specified. If Fs is provided, MAGUNITS must
%   be provided after Fs in the input argument list.
%   
%   The full list of possible values for SPECSTRING (not case sensitive)
%   is:
%
%       'Fp,Fst,Ap,Ast' (minimum order; default)
%       'N,F3dB' 
%       'Nb,Na,F3dB' 
%       'N,F3dB,Fst' (*)  
%       'N,F3dB,Ap' (*)   
%       'N,F3dB,Ast' (*) 
%       'N,F3dB,Ap,Ast' (*)
%       'N,Fc'        
%       'N,Fc,Ap,Ast'   
%       'N,Fp,Ap'     
%       'N,Fp,Ap,Ast'   
%       'N,Fp,F3dB' (*)   
%       'N,Fp,Fst'   
%       'N,Fp,Fst,Ap' (*)  
%       'N,Fp,Fst,Ast' (*) 
%       'N,Fst,Ast'    
%       'N,Fst,Ap,Ast' (*) 
%       'Nb,Na,Fp,Fst' (*)  
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
%   D = FDESIGN.LOWPASS(Fpass, Fstop, Apass, Astop) uses the  default
%   SPECSTRING ('Fp,Fst,Ap,Ast') and sets the passband-edge frequency,
%   stopband-edge frequency, passband ripple, and stopband attenuation.
%
%   % Example #1 - Design a minimum-order equiripple lowpass filter. 
%   % Specify ripples in linear units.
%   d  = fdesign.lowpass('Fp,Fst,Ap,Ast',0.1,0.12,0.1,1e-4,'linear');
%   Hd = design(d, 'equiripple');
%   fvtool(Hd)
%
%   % Example #2 - Design a minimum-phase FIR lowpass filter with
%   % equiripple passband and a stopband that decays as (1/f)^2. (*)
%   d  = fdesign.lowpass('N,Fp,Fst',40,0.3,0.35);
%   designopts(d,'equiripple') % List equiripple design options
%   Hd = design(d,'equiripple','StopbandShape','1/f','StopbandDecay',2,...
%        'minphase',true);
%   fvtool(Hd)
%
%   % Example #3 - Design a Chebyshev Type I IIR filter with a passband 
%   % ripple of 0.5 dB and a 3 dB cutoff frequency at 9600 Hz. (*)
%   Fs = 48000; % Sampling frequency of input signal
%   d  = fdesign.lowpass('N,F3dB,Ap', 10, 9600, .5, Fs);
%   Hd = design(d, 'cheby1');
%   fvtool(Hd)
%
%   %(*) DSP System Toolbox required
%
%   See also FDESIGN, FDESIGN/SETSPECS, FDESIGN/DESIGN, FDESIGN/DESIGNOPTS.

%   Copyright 2004-2013 The MathWorks, Inc.

this = fdesign.lowpass;

[varargin,flag] = finddesignfiltflag(this,varargin);

set(this, 'Response', 'Lowpass');

if flag 
  specObj = this.getcurrentspecs;
  specObj.FromDesignfilt = true;
end

this.setspecs(varargin{:});

capture(this);

% [EOF]
