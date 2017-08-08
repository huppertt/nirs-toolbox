function this = halfband(varargin)
%HALFBAND   Construct a halfband filter designer.
%   D = FDESIGN.HALFBAND constructs a halfband filter designer D.
%
%   D = FDESIGN.HALFBAND('TYPE',TYPE) initializes the filter designer
%   'Type' property with TYPE.  TYPE must be either 'Lowpass' or 'Highpass'
%   and is not case sensitive.
%
%   D = FDESIGN.HALFBAND(SPEC) initializes the filter designer
%   'Specification' property to SPEC.  SPEC is one of the following
%   strings and is not case sensitive:
%
%       'TW,Ast' - (minimum order, default)
%       'N'      
%       'N,Ast' 
%       'N,TW'   
%
%  where 
%       Ast   - Stopband Attenuation (dB)
%       N     - Filter Order
%       TW    - Transition Width
%
%   By default, all frequency specifications are assumed to be in
%   normalized frequency units. Moreover, all magnitude specifications are
%   assumed to be in dB.
%
%   Different specification types may have different design methods
%   available. Use DESIGNMETHODS(D) to get a list of design methods
%   available for a given SPEC.
%
%   D = FDESIGN.HALFBAND(SPEC, SPEC1, SPEC2, ...) initializes the filter
%   designer specifications with SPEC1, SPEC2, etc...
%   Use GET(D, 'DESCRIPTION') for a description of SPEC1, SPEC2, etc.
%
%   D = FDESIGN.HALFBAND(TransitionWidth, Astop) uses the  default
%   SPEC ('TW,Ast') and sets the transition width and stopband        
%   attenuation.
%
%   D = FDESIGN.HALFBAND(...,Fs) specifies the sampling frequency (in
%   Hz). In this case, the transition width, if specified, is also in Hz.
%
%   D = FDESIGN.HALFBAND(...,MAGUNITS) specifies the units for any
%   magnitude specification given in the constructor. MAGUNITS can be one
%   of the following: 'linear', 'dB', or 'squared'. If this argument is
%   omitted, 'dB' is assumed. Note that the magnitude specifications are
%   always converted and stored in dB regardless of how they were
%   specified.
%
%   % Example #1 
%       %Design a minimum order equiripple halfband lowpass filter.
%       %Compare to an elliptic IIR halfband filter
%       d = fdesign.halfband('Type','Lowpass',.01, 80);
%       Hfir = design(d,'equiripple');
%       Hiir = design(d,'ellip');
%       fvtool(Hfir,Hiir)
%
%   % Example #2 
%       %Design a 80th order equiripple halfband highpass filter 
%       %with 70dB of stopband attenuation.
%       d = fdesign.halfband('Type','Highpass','N,Ast',80,70);
%       Hd = design(d,'equiripple');
%
%   % Example #3 
%       %Design a 42th order equiripple halfband lowpass filter 
%       %with a controlled transition width.
%       d = fdesign.halfband('Type','Lowpass','N,TW', 42, .04);
%       designmethods(d);
%       Hd = design(d,'firls');
%
%   For more information about halfband filters, see the
%   <a href="matlab:web([matlabroot,'\toolbox\dsp\dspdemos\html\firhalfbanddemo.html'])">FIR Halfband Filter Design</a> demo. 
%
%   See also FDESIGN, FDESIGN/SETSPECS, FDESIGN/DESIGN.


%   Author(s): J. Schickler
%   Copyright 1999-2010 The MathWorks, Inc.

this = fdesign.halfband;
set(this, 'Response', 'Halfband');
this.setspecs(varargin{:});
capture(this);

% [EOF]
