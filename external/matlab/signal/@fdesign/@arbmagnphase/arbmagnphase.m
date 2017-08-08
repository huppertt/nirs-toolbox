function this = arbmagnphase(varargin)
%ARBMAGNPHASE   Arbitrary Magnitude and Phase filter designer.
%   D = FDESIGN.ARBMAGNPHASE constructs an arbitrary magnitude filter designer D.
%
%   D = FDESIGN.ARBMAGNPHASE(SPEC) initializes the filter designer
%   'Specification' property to SPEC.  SPEC is one of the following
%   strings and is not case sensitive:
%
%       'N,F,H'       - Single-band design (default)
%       'Nb,Na,F,H'   - Single-band IIR design
%       'N,B,F,H'     - Multi-band design
%
%  where 
%       H  - Complex Frequency Response 
%       B  - Number of Bands
%       F  - Frequency Vector
%       N  - Filter Order
%       Nb - Numerator Order
%       Na - Denominator Order
%
%   By default, all frequency specifications are assumed to be in
%   normalized frequency units. 
%
%   Different specification types may have different design methods
%   available. Use DESIGNMETHODS(D) to get a list of design methods
%   available for a given SPEC.
%
%   D = FDESIGN.ARBMAGNPHASE(SPEC, SPEC1, SPEC2, ...) initializes the filter
%   designer specifications with SPEC1, SPEC2, etc. 
%   Use GET(D, 'DESCRIPTION') for a description of SPEC1, SPEC2, etc.
%
%   D = FDESIGN.ARBMAGNPHASE(N, F, H) uses the  default SPEC ('N,F,H') and
%   sets the order, the frequency vector, and the frequency response vector.
%
%   D = FDESIGN.ARBMAGNPHASE(...,Fs) specifies the sampling frequency (in Hz).
%   In this case, all other frequency specifications are also in Hz.
%
%   % Example #1 - Model a Complex Analog Filter.
%   d=fdesign.arbmagnphase('N,F,H',100);
%   design(d,'freqsamp');
%
%   % Example #2 - Design a Bandpass Filter with a Low Group Delay
%   N  = 50;  % Group Delay of linear phase filter would be 25
%   gd = 12; % Desired Group Delay
%   F1 = linspace(0,.25,30); F2=linspace(.3,.56,40); F3=linspace(.62,1,30);
%   H1 = zeros(size(F1)); H2 = exp(-1i*pi*gd*F2); H3 = zeros(size(F3));
%   d  = fdesign.arbmagnphase('N,B,F,H',50,3,F1,H1,F2,H2,F3,H3); 
%   Hd = design(d,'equiripple');
%   fvtool(Hd)
%
%   For more information, see the <a href="matlab:web([matlabroot,'\toolbox\dsp\dspdemos\html\arbmagnphasedemo.html'])">Arbitrary Magnitude and Phase Demo</a>. 
%
%   See also FDESIGN, FDESIGN/SETSPECS, FDESIGN/DESIGN.


%   Author(s): V. Pellissier
%   Copyright 2005-2010 The MathWorks, Inc.

this = fdesign.arbmagnphase;

set(this, 'Response', 'Arbitrary Magnitude and Phase');

this.setspecs(varargin{:});

capture(this);

% [EOF]
