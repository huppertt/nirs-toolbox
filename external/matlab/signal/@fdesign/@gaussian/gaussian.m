function this = gaussian(varargin)
%GAUSSIAN   Construct a Gaussian pulse shaping filter designer.
%   D = FDESIGN.GAUSSIAN(SPS,SPECSTRING,VALUE1,VALUE2,...)
%   constructs a Gaussian pulse shaping filter designer D. Note that D is not
%   the design itself, it only contains the design specifications. In order to
%   design the filter, one needs to invoke the DESIGN method on D. 
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
%   D = FDESIGN.GAUSSIAN(...,Fs) provides the sampling frequency of the
%   signal to be filtered. Fs must be specified as a scalar trailing the
%   other numerical values provided. For this case, Fs is assumed to be in
%   Hz and used for analysis and visualization purposes.
%
%   D = FDESIGN.GAUSSIAN(...,MAGUNITS) specifies the units for any magnitude
%   specification given. MAGUNITS can be one of the following: 'linear', 'dB',
%   or 'squared'. If this argument is omitted, 'dB' is assumed. Note that the
%   magnitude specifications are always converted and stored in dB regardless of
%   how they were specified. If Fs is provided, MAGUNITS must be provided after
%   Fs in the input argument list.
%   
%   The full list of possible values for SPECSTRING (not case sensitive)
%   is:
%           'Nsym,BT' 
%
%  where 
%       BT    - 3 dB Bandwidth - Symbol Time product
%       Nsym  - Filter Order in symbols
%
%   D = FDESIGN.GAUSSIAN(sps, Nsym, BT) uses the  default
%   SPECSTRING.
%
%   % Example - Design a gaussian pulse shaping FIR filter of order 8 symbols, 
%   %  bandwidth-symbol time product of 0.50, and 10 samples per symbol.
%   h  = fdesign.gaussian(10,'Nsym,BT',8,0.50);
%   Hd = design(h);
%   fvtool(Hd)
%  
%   See also FDESIGN, FDESIGN/SETSPECS, FDESIGN/DESIGN, FDESIGN/DESIGNOPTS.

%   Copyright 2008 The MathWorks, Inc.

this = fdesign.gaussian;

set(this, 'Response', 'Gaussian');

this.setspecs(varargin{:});

capture(this);

% [EOF]
