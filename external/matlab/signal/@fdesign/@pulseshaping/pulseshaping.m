function this = pulseshaping(varargin)
%PULSESHAPING   Construct a pulse shaping filter designer.
%
%   WARNING: fdesign.pulseshaping is not recommended. Use rcosdesign or
%            gaussdesign instead. 
%
%   D = FDESIGN.PULSESHAPING(SPS,SHAPE,SPECSTRING,VALUE1,VALUE2,...)
%   constructs a pulse shaping filter designer D. Note that D is not the design
%   itself, it only contains the design specifications. In order to design the
%   filter, one needs to invoke the DESIGN method on D.
%   For example (more examples below):
%   D = fdesign.pulseshaping(8,'Raised Cosine','Nsym,Beta',6,0.25);
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
%   D = FDESIGN.PULSESHAPING(...,Fs) provides the sampling frequency of the
%   signal to be filtered. Fs must be specified as a scalar trailing the
%   other numerical values provided. For this case, Fs is assumed to be in
%   Hz and used for analysis and visualization purposes.
%
%   D = FDESIGN.PULSESHAPING(...,MAGUNITS) specifies the units for any magnitude
%   specification given. MAGUNITS can be one of the following: 'linear', 'dB',
%   or 'squared'. If this argument is omitted, 'dB' is assumed. Note that the
%   magnitude specifications are always converted and stored in dB regardless of
%   how they were specified. If Fs is provided, MAGUNITS must be provided after
%   Fs in the input argument list.
%
%   The full list of possible values for SPECSTRING (not case sensitive)
%   is:
%       For PULSESHAPE 'Raised Cosine' and 'Square Root Raised Cosine':
%           'Ast,Beta' (minimum order; default)
%           'Nsym,Beta'
%           'N,Beta'
%       For PULSESHAPE 'Gaussian':
%           'Nsym,BT' (default)
%
%  where
%       Ast   - Stopband Attenuation (dB)
%       Beta  - Rolloff factor
%       Nsym  - Filter Order in symbols (must be even for raised cosine filters)
%       N     - Filter Order (must be even)
%       BT    - Bandwidth - Symbol Time product
%
%   D = FDESIGN.PULSESHAPING(sps, shape, Astop, Beta) uses the  default
%   SPECSTRING.
%
%   % Example #1 - Design a raised cosine windowed FIR filter with stop band
%   % attenuation of 60dB, rolloff factor of 0.50, and 8 samples
%   % per symbol.
%   h  = fdesign.pulseshaping(8,'Raised Cosine','Ast,Beta',60,0.50);
%   Hd = design(h);
%   fvtool(Hd)
%
%   % Example #2 - Design a raised cosine windowed FIR filter of order 8 symbols,
%   %  rolloff factor of 0.50, and 10 samples per symbol.
%   h  = fdesign.pulseshaping(10,'Raised Cosine','Nsym,Beta',8,0.50);
%   Hd = design(h);
%   fvtool(Hd)
%
%   % Example #3 - Design a square root raised cosine windowed FIR filter of order
%   % 42, rolloff factor of 0.25, and 10 samples per symbol.
%   h  = fdesign.pulseshaping(10,'Square Root Raised Cosine','N,Beta',42);
%   Hd = design(h);
%   fvtool(Hd)
%
%   % Example #4 - Design a Gaussian windowed FIR filter of order 3 symbols, 
%   % bandwidth-symbol time product of 0.4, and 10 samples per symbol.
%   h  = fdesign.pulseshaping(10,'Gaussian','Nsym,BT',3,0.4);
%   Hd = design(h);
%   fvtool(Hd)
%
%
%   See also FDESIGN, FDESIGN/SETSPECS, FDESIGN/DESIGN, FDESIGN/DESIGNOPTS.

%   Copyright 2008-2011 The MathWorks, Inc.


this = fdesign.pulseshaping;

if nargin > 0 && ~isnumeric(varargin{1})
  error(message('signal:fdesign:abstractpulseshape:setspecs:invalidInputSingleRate'));
end

if nargin < 2
    this.PulseShape = 'Raised Cosine';
else
    this.PulseShape = varargin{2};
    varargin(2) = [];
end

setspecs(this.PulseShapeObj, varargin{:})

% [EOF]
