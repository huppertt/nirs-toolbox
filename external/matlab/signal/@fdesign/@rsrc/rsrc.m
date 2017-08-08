function this = rsrc(L, M, DT, varargin)
%RSRC Construct a rational sample-rate converter (rsrc) filter designer.
%   D = FDESIGN.RSRC(L, M) constructs an rsrc filter designer D with
%   an 'InterpolationFactor' of L and a 'DecimationFactor' of M.  If L is
%   not specified, it defaults to 3.  If M is not specified it defaults to
%   2.
%
%   D = FDESIGN.RSRC(L, M, , RESPONSE) initializes the filter designer
%   'Response' property with RESPONSE.  RESPONSE is one of the following
%   strings and is not case sensitive:
%       'Nyquist' (default)
%       'Halfband'
%       'Lowpass'
%       'CIC Compensator'
%       'Inverse-sinc Lowpass'
%       'Inverse-sinc Highpass'
%       'Highpass'
%       'Hilbert'
%       'Bandpass'
%       'Bandstop'
%       'Differentiator'
%       'Arbitrary Magnitude'
%       'Arbitrary Magnitude and Phase'
%       'Raised Cosine'
%       'Square Root Raised Cosine'
%       'Gaussian'
%
%   D = FDESIGN.RSRC(L, M, RESPONSE, SPEC) initializes the filter
%   designer 'Specification' property with SPEC.  Use SET(D, 'SPECIFICATION')
%   to get all available specifications for the response RESPONSE.
%
%   Different design and specification types will have different design
%   methods available. Use DESIGNMETHODS(D) to get a list of design
%   methods available for a given SPEC.
%
%   D = FDESIGN.RSRC(L, M, RESPONSE, SPEC, SPEC1, SPEC2, ...) initializes
%   the filter designer specifications with SPEC1, SPEC2, etc.
%   Use GET(D, 'DESCRIPTION') for a description of SPEC1, SPEC2, etc.
%
%   D = FDESIGN.RSRC(...,Fs) specifies the sampling frequency (in Hz). In
%   this case, all frequency properties are also in Hz.
%
%   D = FDESIGN.RSRC(...,MAGUNITS) specifies the units for any magnitude
%   specification given in the constructor. MAGUNITS can be one of the
%   following: 'linear', 'dB', or 'squared'. If this argument is omitted,
%   'dB' is assumed. Note that the magnitude specifications are always
%   converted and stored in dB regardless of how they were specified.
%
%   Some responses have additional constructor inputs. The syntaxes for
%   these special cases are listed below:
%
%   D = FDESIGN.RSRC(L, M, 'CIC compensator', DD, N, R, SPEC, ...) creates
%   a CIC compensator rational sample-rate converter filter designer with
%   the 'InterpolationFactor' property set to L, the 'DecimationFactor'
%   property set to M, the 'DifferentialDelay' property set to DD, the
%   'NumberOfSections' property set to N, and the 'CICRateChangeFactor'
%   property set to R. In general, the design method implements a filter
%   with a passband response equal to an inverse Dirichlet sinc that
%   matches exactly the inverse passband response of a CIC filter with a
%   differential delay equal to DD, number of sections equal to N, and rate
%   change factor equal to R. When R is not specified or is set to 1, the
%   design method implements a filter with a passband response equal to an
%   inverse sinc that is an approximation to the true inverse passband
%   response of the CIC filter.
%
%   D = FDESIGN.RSRC(L, M, PULSESHAPERESP, SPS, SPEC, ...) where
%   PULSESHAPERESP equals 'Raised Cosine', 'Square Root Raised Cosine', or
%   'Gaussian', creates a pulse shaping rational sample-rate converter
%   filter designer with the 'InterpolationFactor' property set to L, the
%   'DecimationFactor' property set to M, and the 'SamplesPerSymbol'
%   property set to SPS.
%
%   % Example #1 - Design a minimum order Nyquist sample-rate converter.
%   d = fdesign.rsrc(5, 3, 'Nyquist',5,0.05, 40);
%   designmethods(d)
%   hm = design(d, 'kaiserwin');
%
%   % Example #2 - Design a 30th order Nyquist sample-rate converter.
%   d = fdesign.rsrc(5, 3, 'Nyquist', 5, 'N,TW', 30)
%   design(d)
%
%   % Example #3 - Specify frequencies in Hz.
%   d = fdesign.rsrc(5, 3, 'Nyquist', 5, 'N,TW', 12, 0.1, 5)
%   designmethods(d);
%   design(d, 'equiripple');
%
%   % Example #4 - Specify a stopband ripple in linear units
%   d = fdesign.rsrc(4,7,'Nyquist',7,'TW,Ast',.1,1e-3,5,'linear') % 1e-3 = 60dB
%   design(d)
%
%   See also FDESIGN, FDESIGN/SETSPECS, FDESIGN/DESIGN.

%   Copyright 2005-2011 The MathWorks, Inc.

this = fdesign.rsrc;

set(this, 'MultirateType', 'Rational Sample Rate Converter');

needsDefaults = true;

if nargin > 0
    set(this, 'InterpolationFactor', L);
    if nargin > 1
        set(this, 'DecimationFactor', M);
        if nargin > 2
            if isa(DT, 'fdesign.abstracttypewspecs')
                set(this, 'AllFDesign', copy(DT), 'CurrentFDesign', []);
                [~, DT] = strtok(class(DT), '.');
                DT(1) = [];
                needsDefaults = false;
            end
            newresp = mapresponse(this, DT);
            if strcmpi(newresp, this.Response)
                updatecurrentfdesign(this);
            else
                set(this, 'Response', newresp);
            end
        end
    end
end

if needsDefaults
    multiratedefaults(this.CurrentFDesign, max(getratechangefactors(this)));
end

setspecs(this, varargin{:});

% [EOF]
