function this = nyquist(varargin)
%NYQUIST   Construct a nyquist filter designer.
%   D = FDESIGN.NYQUIST(L) constructs a Lth band Nyquist filter designer D.
%
%   The band of a Nyquist filter is the inverse of the cutoff frequency in
%   terms of normalized units. For instance, a 4th-band filter has a cutoff
%   of 1/4. The case L=2 is referred to as a halfband filter. Use
%   FDESIGN.HALFBAND for more options with halfband filter design.
%
%   D = FDESIGN.NYQUIST(L,SPEC) initializes the filter designer
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
%   D = FDESIGN.NYQUIST(L, SPEC, SPEC1, SPEC2, ...) initializes the filter
%   designer specifications with SPEC1, SPEC2, etc...
%   Use GET(D, 'DESCRIPTION') for a description of SPEC1, SPEC2, etc.
%
%   D = FDESIGN.NYQUIST(L, TransitionWidth, Astop) uses the  default
%   SPEC ('TW,Ast') and sets the transition width and stopband        
%   attenuation.
%
%   D = FDESIGN.NYQUIST(...,Fs) specifies the sampling frequency (in
%   Hz). In this case, the transition width, if specified, is also in Hz.
%
%   D = FDESIGN.NYQUIST(...,MAGUNITS) specifies the units for any
%   magnitude specification given in the constructor. MAGUNITS can be one
%   of the following: 'linear', 'dB', or 'squared'. If this argument is
%   omitted, 'dB' is assumed. Note that the magnitude specifications are
%   always converted and stored in dB regardless of how they were
%   specified.
%
%   % Example #1 - Design a minimum order, 4th band Nyquist filter.
%   d = fdesign.nyquist(4,'TW,Ast',.01, 80);
%   designmethods(d);
%   Hd = design(d, 'kaiserwin'); 
%   fvtool(Hd)
%
%   % Example #2 - Design a 42nd order, 5th band Nyquist filter.
%   d = fdesign.nyquist(5,'N,Ast', 42, 80)
%   design(d)
%
%   % Example #3 - Control the transition width.
%   d = fdesign.nyquist(5,'N,TW', 42, .1)
%   design(d)
%
%   % Example #4 - Design a 2nd band (halfband) Nyquist filter. Compare FIR
%   % equiripple and IIR elliptic designs
%   d = fdesign.nyquist(2,'TW,Ast',.1,80);
%   H(1) = design(d,'equiripple');
%   H(2) = design(d,'ellip');
%   hfvt = fvtool(H); legend(hfvt,'Equiripple','Elliptic')
%
%   See also FDESIGN, FDESIGN/SETSPECS, FDESIGN/DESIGN.

%   Author(s): J. Schickler
%   Copyright 1999-2005 The MathWorks, Inc.

this = fdesign.nyquist;

set(this, 'Response', 'Nyquist');

if nargin>0,
    defaulttw(this,varargin{1})
end

this.setspecs(varargin{:});
capture(this);

% [EOF]
