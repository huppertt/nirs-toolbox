function this = hilbert(varargin)
%HILBERT   Construct a Hilbert filter designer.
%   D = FDESIGN.HILBERT constructs a Hilbert filter designer D.
%
%   D = FDESIGN.HILBERT(SPEC) initializes the filter designer
%   'Specification' property to SPEC.  SPEC is one of the following
%   strings and is not case sensitive:
%
%       'N,TW'  - (default)
%       'TW,Ap' - Minimum order (*)
%
%   Different specification may have different design methods available.
%   Use DESIGNMETHODS(D) to get a list of design methods available for a
%   given SPEC.
%
%   D = FDESIGN.HILBERT(SPEC, SPEC1, SPEC2, ...) initializes the filter
%   designer specifications with SPEC1, SPEC2, etc...
%   Use GET(D, 'DESCRIPTION') for a description of SPEC1, SPEC2, etc.
%
%   D = FDESIGN.HILBERT(...,Fs) specifies the sampling frequency
%   (in Hz). In this case, all other frequency specifications are also in
%   Hz.
%
%   D = FDESIGN.HILBERT(...,MAGUNITS) specifies the units for any
%   magnitude specification given in the constructor. MAGUNITS can be one
%   of the following: 'linear', 'dB', or 'squared'. If this argument is
%   omitted, 'dB' is assumed. Note that the magnitude specifications are
%   always converted and stored in dB regardless of how they were
%   specified.
%
%   % Example #1 - Design a 30th order type III Hilbert Transformer.
%   d = fdesign.hilbert(30,.2);
%   designmethods(d);
%   Hd = design(d,'firls'); 
%   fvtool(Hd,'MagnitudeDisplay','Zero-phase','FrequencyRange','[-pi, pi)')
%
%   % Example #2 - Design a 35th order type IV Hilbert Transformer.
%   d = fdesign.hilbert('N,TW',35,.1);
%   Hd = design(d,'equiripple'); 
%   fvtool(Hd,'MagnitudeDisplay','Zero-phase','FrequencyRange','[-pi, pi)')
%
%   % Example #3 - Design a minimum-order Hilbert Transformer with a
%   % sampling frequency of 100, compare two IIR designs. (*)
%   d = fdesign.hilbert('TW,Ap',1,.1,100);
%   H(1) = design(d,'iirlinphase');
%   H(2) = design(d,'ellip');
%   hfvt = fvtool(H); legend(hfvt,'Linear phase IIR','Elliptic IIR')
%
%   %(*) DSP System Toolbox required
%
%   See also FDESIGN, FDESIGN/SETSPECS, FDESIGN/DESIGN.


%   Copyright 2005-2013 The MathWorks, Inc.

this = fdesign.hilbert;

[varargin,flag] = finddesignfiltflag(this,varargin);

set(this, 'Response', 'Hilbert Transformer');

if flag 
  specObj = this.getcurrentspecs;
  specObj.FromDesignfilt = true;
end

this.setspecs(varargin{:});

capture(this);


% [EOF]
