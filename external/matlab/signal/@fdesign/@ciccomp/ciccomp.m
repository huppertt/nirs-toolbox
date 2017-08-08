function this = ciccomp(M, N, varargin)
%CICCOMP   Construct a CIC compensator filter designer.
%   D = FDESIGN.CICCOMP constructs a CIC compensator filter designer D.
%
%   D = FDESIGN.CICCOMP(DELAY, NSECTIONS, RCIC) constructs a CIC
%   compensator filter designer with DifferentialDelay set to DELAY,
%   NumberOfSections set to NSECTIONS, and CICRateChangeFactor set to RCIC.
%   By default these parameters are equal to 2, 1, and 1 respectively. In
%   general, the design method implements a filter with a passband response
%   equal to an inverse Dirichlet sinc that matches exactly the inverse
%   passband response of a CIC filter with a differential delay equal to
%   DELAY, number of sections equal to NSECTIONS, and rate change factor
%   equal to RCIC. When RCIC is not specified or is set to 1, the design
%   method implements a filter with a passband response equal to an inverse
%   sinc that is an approximation to the true inverse passband response of
%   the CIC filter.
%
%   D = FDESIGN.CICCOMP(..., SPEC) constructs a CIC compensator and
%   sets its 'Specification' property to SPEC.  SPEC is not case
%   sensitive and must be one of the following:
%
%       'Fp,Fst,Ap,Ast' (minimum order, default)
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
%   D = FDESIGN.CICCOMP(..., SPEC, SPEC1, SPEC2, ...) initializes the
%   filter designer specifications in the order they are specified in the
%   SPEC input. Use GET(D, 'DESCRIPTION') for a description of SPEC1,
%   SPEC2, etc.
%
%   D = FDESIGN.CICCOMP(...,Fs) provides the sampling frequency of the
%   signal to be filtered. Fs must be specified as a scalar trailing the
%   other numerical values provided. For this case, Fs is assumed to be in
%   Hz as are all other frequency values provided. Note that you don't
%   change the specification string in this case. 
%
%   D = FDESIGN.CICCOMP(...,MAGUNITS) specifies the units for any magnitude
%   specification given. MAGUNITS can be one of the following: 'linear',
%   'dB', or 'squared'. If this argument is omitted, 'dB' is assumed. Note
%   that the magnitude specifications are always converted and stored in dB
%   regardless of how they were specified. If Fs is provided, MAGUNITS must
%   be provided after Fs in the input argument list.
%
%   % Example #1 - Design a minimum-order CIC compensator that compensates
%   %              for the droop in the passband of the CIC decimator.
%   Fs = 96e3;   % Input sampling frequency 
%   Fpass = 4e3; % Frequency band of interest
%   M = 6;       % Decimation factor of CIC      
%   Hcic = design(fdesign.decimator(M,'CIC',1,Fpass,60,Fs));
%   Hd = cascade(dfilt.scalar(1/gain(Hcic)),Hcic);
%   d = fdesign.ciccomp(Hcic.DifferentialDelay, Hcic.NumberOfSections,...
%         Hcic.DecimationFactor,Fpass,4.5e3,.1,60,Fs/M);
%   Hd(2) = design(d);
%   hfvt=fvtool(Hd(1),Hd(2),cascade(Hd(1),Hd(2)),'Fs',[96e3 96e3/M 96e3], ...
%         'ShowReference', 'off');
%   legend(hfvt, 'CIC decimator','CIC compensator','Overall response', ...
%         'Location', 'Northeast');
%
%   % Example #2 - Design a compensator whose stopband decays like (1/f)^2.
%   Hd(3) = design(d, 'equiripple', 'StopbandShape', '1/f', 'StopbandDecay', 2);
%   hfvt=fvtool(Hd(1),Hd(3),cascade(Hd(1),Hd(3)),'Fs',[96e3 96e3/M 96e3], ...
%         'ShowReference', 'off');
%   legend(hfvt, 'CIC decimator','CIC compensator','Overall response', ...
%         'Location', 'Northeast');
%
%   See also FDESIGN, FDESIGN/SETSPECS, FDESIGN/DESIGN.

%   Copyright 2005-2011 The MathWorks, Inc.

this = fdesign.ciccomp;

if nargin < 1, M = 1; end
if nargin < 2, N = 2; end

set(this, 'Response', 'CIC Compensator');

this.setspecs(M, N, varargin{:});

capture(this);

% [EOF]
