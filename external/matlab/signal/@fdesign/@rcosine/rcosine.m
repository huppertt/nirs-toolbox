function this = rcosine(varargin)
%RCOSINE   Construct a raised cosine pulse shaping filter designer.
%   D = FDESIGN.RCOSINE(SPS,SPECSTRING,VALUE1,VALUE2,...)
%   constructs a raised cosine pulse shaping filter designer D. Note that D is
%   not the design itself, it only contains the design specifications. In order
%   to design the filter, one needs to invoke the DESIGN method on D. 
%   For example (more examples below): 
%   D = fdesign.rcosine(8,'Nsym,Beta',6,0.25); 
%   H = design(D); % H is a DFILT
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
%   corresponding specification. In the example above, this means that Nsym =
%   6, Beta (Rolloff factor) = 0.25. Use get(D, 'description') for a
%   description of VALUE1, VALUE2, etc.
%
%   By default, all frequency specifications are assumed to be in
%   normalized frequency units. Moreover, all magnitude specifications are
%   assumed to be in dB.
%
%   D = FDESIGN.RCOSINE(...,Fs) provides the sampling frequency of the
%   signal to be filtered. Fs must be specified as a scalar trailing the
%   other numerical values provided. For this case, Fs is assumed to be in
%   Hz and used for analysis and visualization purposes.
%
%   D = FDESIGN.RCOSINE(...,MAGUNITS) specifies the units for any magnitude
%   specification given. MAGUNITS can be one of the following: 'linear', 'dB',
%   or 'squared'. If this argument is omitted, 'dB' is assumed. Note that the
%   magnitude specifications are always converted and stored in dB regardless of
%   how they were specified. If Fs is provided, MAGUNITS must be provided after
%   Fs in the input argument list.
%   
%   The full list of possible values for SPECSTRING (not case sensitive)
%   is:
%           'Ast,Beta' (minimum order; default)
%           'Nsym,Beta' 
%           'N,Beta'
%
%  where 
%       Ast   - Stopband Attenuation (dB)
%       Beta  - Rolloff factor
%       Nsym  - Filter Order in symbols (must be even)
%       N     - Filter Order (must be even)
%
%   D = FDESIGN.RCOSINE(sps, Astop, Beta) uses the  default
%   SPECSTRING.
%
%   % Example #1 - Design a raised cosine windowed FIR filter with stop band 
%   % attenuation of 60dB, rolloff factor of 0.50, and 8 samples 
%   % per symbol.
%   h  = fdesign.rcosine(8,'Ast,Beta',60,0.50);
%   Hd = design(h);
%   fvtool(Hd)
%
%   % Example #2 - Design a raised cosine windowed FIR filter of order 8 symbols, 
%   %  rolloff factor of 0.50, and 10 samples per symbol.
%   h  = fdesign.rcosine(10,'Nsym,Beta',8,0.50);
%   Hd = design(h);
%   fvtool(Hd)
%  
%   See also FDESIGN, FDESIGN/SETSPECS, FDESIGN/DESIGN, FDESIGN/DESIGNOPTS.

%   Copyright 2008 The MathWorks, Inc.

this = fdesign.rcosine;

set(this, 'Response', 'Raised Cosine');

this.setspecs(varargin{:});

capture(this);

% [EOF]
