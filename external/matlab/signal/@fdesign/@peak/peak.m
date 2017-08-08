function this = peak(varargin)
%PEAK    Construct a peaking filter designer.
%   D = FDESIGN.PEAK(SPECSTRING,VALUE1,VALUE2,...) constructs a
%   peaking filter designer D. Note that D is not the design
%   itself, it only contains the design specifications. In order to design
%   the filter, one needs to invoke the DESIGN method on D.
%   For example (more examples below):
%   D = fdesign.peak('N,F0,Q',4,0.3,6.5);
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
%   corresponding specification. In the example above, this means that N =
%   4, F0 = 0.3, and Q = 6.5. Use get(D,'description') for a description of
%   VALUE1, VALUE2, etc.
%
%   By default, all frequency specifications are assumed to be in
%   normalized frequency units. Moreover, all magnitude specifications are
%   assumed to be in dB.
%
%   D = FDESIGN.PEAK(...,Fs) provides the sampling frequency of the
%   signal to be filtered. Fs must be specified as a scalar trailing the
%   other numerical values provided. For this case, Fs is assumed to be in
%   Hz as are all other frequency values provided. Note that you don't
%   change the specification string in this case. In the example above, if
%   the input signal is sampled at 8 kHz, we can obtain the same filter by
%   specifying the frequencies in Hz as:
%   D = fdesign.peak('N,F0,Q',4,1200,6.5,8e3);
%   H = design(D);
%
%   D = FDESIGN.PEAK(...,MAGUNITS) specifies the units for any magnitude
%   specification given. MAGUNITS can be one of the following: 'linear',
%   'dB', or 'squared'. If this argument is omitted, 'dB' is assumed. Note
%   that the magnitude specifications are always converted and stored in dB
%   regardless of how they were specified. If Fs is provided, MAGUNITS must
%   be provided after Fs in the input argument list.
%
%   The full list of possible values for SPECSTRING (not case sensitive) is:
%
%      'N,F0,Q' (default)
%      'N,F0,Q,Ap'
%      'N,F0,Q,Ast'
%      'N,F0,Q,Ap,Ast'
%      'N,F0,BW'
%      'N,F0,BW,Ap'
%      'N,F0,BW,Ast'
%      'N,F0,BW,Ap,Ast'
%
%   where 
%       N   - Filter Order (must be even)
%       F0  - Center Frequency
%       Q   - Quality Factor
%       BW  - 3-dB Bandwidth
%       Ap  - Passband Ripple (dB)
%       Ast - Stopband Attenuation (dB)
%
%   D = FDESIGN.PEAK(n,f0,q) uses the  default SPECSTRING ('N,F0,Q') and
%   sets the corresponding values to n, f0, and q.
%
%   % Example #1 - Design a peaking filter with a passband ripple of 1 dB.
%   d  = fdesign.peak('N,F0,Q,Ap',6,0.5,10,1);
%   Hd = design(d);
%   fvtool(Hd)
%
%   % Example #2 - Design a Chebyshev Type II peaking filter with a
%   % stopband attenuation of 80 dB
%   d = fdesign.peak('N,F0,BW,Ast',8,.65,.02,80);
%   Hd = design(d,'cheby2');
%   fvtool(Hd)
%
%   See also FDESIGN, FDESIGN/SETSPECS, FDESIGN/DESIGN, FDESIGN/DESIGNOPTS.


%   Author(s): S. Orfanidis and R. Losada
%   Copyright 2006 The MathWorks, Inc.

this = fdesign.peak;

set(this, 'Response', 'Peaking Filter');

this.setspecs(varargin{:});

capture(this);

% [EOF]
