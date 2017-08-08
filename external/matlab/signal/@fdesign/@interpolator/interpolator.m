function this = interpolator(L, DT, varargin)
%INTERPOLATOR   Construct an interpolator filter designer.
%   D = FDESIGN.INTERPOLATOR(L) constructs an interpolator filter designer
%   D with an 'InterpolationFactor' of L. If L is not specified, it
%   defaults to 2.
%
%   D = FDESIGN.INTERPOLATOR(L, RESPONSE) initializes the filter designer
%   'Response' property with RESPONSE.  RESPONSE is one of the following
%   strings and is not case sensitive: 
%       'Nyquist' (default)
%       'Halfband'
%       'Lowpass'
%       'CIC'
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
%   D = FDESIGN.INTERPOLATOR(L, RESPONSE, SPEC) initializes the filter
%   designer 'Specification' property with SPEC.  Use SET(D, 'SPECIFICATION')
%   to get all available specifications for the response RESPONSE. 
%
%   Different design and specification types will have different design
%   methods available. Use DESIGNMETHODS(D) to get a list of design
%   methods available for a given SPEC.
%
%   D = FDESIGN.INTERPOLATOR(L, RESPONSE, SPEC, SPEC1, SPEC2, ...)
%   initializes the filter designer specifications with SPEC1, SPEC2, etc.
%   Use GET(D, 'DESCRIPTION') for a description of SPEC1, SPEC2, etc.
%
%   D = FDESIGN.INTERPOLATOR(...,Fs) specifies the sampling frequency (in
%   Hz). In this case, all frequency properties are also in Hz.
%
%   D = FDESIGN.INTERPOLATOR(...,MAGUNITS) specifies the units for any
%   magnitude specification given in the constructor. MAGUNITS can be one
%   of the following: 'linear', 'dB', or 'squared'. If this argument is
%   omitted, 'dB' is assumed. Note that the magnitude specifications are
%   always converted and stored in dB regardless of how they were
%   specified.
%
%   Some responses have additional constructor inputs. The syntaxes for
%   these special cases are listed below:
%
%   D = FDESIGN.INTERPOLATOR(L, 'CIC', DD, SPEC, ...) creates a CIC
%   interpolator filter designer with the 'InterpolationFactor' property
%   set to L, and the 'DifferentialDelay' property set to DD.
%
%   D = FDESIGN.INTERPOLATOR(L, 'CIC compensator', DD, N, R, SPEC, ...)
%   creates a CIC compensator interpolator filter designer with the
%   'InterpolationFactor' property set to L, the 'DifferentialDelay'
%   property set to DD, the 'NumberOfSections' property set to N, and the
%   'CICRateChangeFactor' property set to R. In general, the design method
%   implements a filter with a passband response equal to an inverse
%   Dirichlet sinc that matches exactly the inverse passband response of a
%   CIC filter with a differential delay equal to DD, number of sections
%   equal to N, and rate change factor equal to R. When R is not specified
%   or is set to 1, the design method implements a filter with a passband
%   response equal to an inverse sinc that is an approximation to the true
%   inverse passband response of the CIC filter.
%
%   D = FDESIGN.INTERPOLATOR(L, PULSESHAPERESP, SPS, SPEC, ...) where
%   PULSESHAPERESP equals 'Raised Cosine', 'Square Root Raised Cosine', or
%   'Gaussian', creates a pulse shaping interpolator filter designer with
%   the 'InterpolationFactor' property set to L, and the 'SamplesPerSymbol'
%   property set to SPS.
%
%   % Example #1 - Design a CIC interpolator for a signal sampled at 19200 Hz
%   % with a differential delay of two and that attenuates images beyond 50 Hz
%   % by at least 80 dB.
%   DD  = 2;     % Differential delay
%   Fp  = 50;    % Passband of interest
%   Ast = 80;    % Minimum attenuation of alias components in passband
%   Fs  = 600;   % Sampling frequency for input signal
%   L   = 32;    % Interpolation factor
%   d   = fdesign.interpolator(L,'cic',DD,'Fp,Ast',Fp,Ast,L*Fs);
%   hm  = design(d);
%
%   % Example #2 - Design a minimum-order CIC compensator that interpolates by
%   % 4 and pre-compensates for the droop in the passband for the CIC from the
%   % previous example. 
%   Nsecs = hm.NumberOfSections;
%   R     = hm.InterpolationFactor;
%   d     = fdesign.interpolator(4,'ciccomp',DD,Nsecs,R,'Fp,Fst,Ap,Ast',50,100,0.1,80,Fs);
%   hmc   = design(d, 'equiripple');
%   hmc.Arithmetic = 'fixed';
%   % Analyze filters individually plus compound response
%   hfvt = fvtool(hmc,hm,cascade(hmc,hm),'Fs',[Fs,L*Fs,L*Fs],'ShowReference','off');
%   legend(hfvt,'CIC pre-compensator','CIC interpolator','Overall response');
%
%   % Example #3 - Design a minimum-order Nyquist interpolator using a Kaiser
%   % window. Compare to a multistage design. 
%   L   = 15;   % Interpolation factor. Also the Nyquist band.
%   TW  = 0.05; % Normalized transition width
%   Ast = 40;   % Minimum stopband attenuation in dB
%   d = fdesign.interpolator(L,'nyquist',L,TW,Ast);
%   hm  = design(d,'kaiserwin');  
%   hm2 = design(d,'multistage');
%   hfvt = fvtool(hm,hm2); legend(hfvt,'Kaiser window','Multistage');
%
%   % Example #4 - Design a lowpass interpolator for an interpolation factor of 8.
%   % Compare a single-stage equiripple design with multistage designs.
%   L = 8; % Interpolation factor;
%   d = fdesign.interpolator(L,'lowpass');
%   hm(1) = design(d,'equiripple');
%   hm(2) = design(d,'multistage','Usehalfbands',true); % Use halfband filters whenever possible
%   hm(3) = design(d,'multistage','Usehalfbands',true,...
%    'HalfbandDesignMethod','iirlinphase'); % Use quasi linear-phase IIR halfbands
%   hfvt = fvtool(hm);
%   legend(hfvt,'Single-stage equiripple','Multistage','Multistage with IIR halfbands')
%
%   See also FDESIGN, FDESIGN/SETSPECS, FDESIGN/DESIGN.

%   Copyright 2005-2011 The MathWorks, Inc.


this = fdesign.interpolator;

set(this, 'MultirateType', 'Interpolator');

needsDefaults = true;

if nargin > 0
    set(this, 'InterpolationFactor', L);
    if nargin > 1
        if isa(DT, 'fdesign.pulseshaping')
            error(message('signal:fdesign:interpolator:interpolator:unsupportedResponse', class( DT )));
        end
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

if needsDefaults
    multiratedefaults(this.CurrentFDesign, this.InterpolationFactor);
end

setspecs(this, varargin{:});

% [EOF]
