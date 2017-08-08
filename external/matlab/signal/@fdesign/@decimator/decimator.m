function this = decimator(M, DT, varargin)
%DECIMATOR   Construct a decimator filter designer.
%   D = FDESIGN.DECIMATOR(M) constructs a decimator filter designer D with
%   a 'DecimationFactor' of M.  If M is not specified, it defaults to 2.
%
%   D = FDESIGN.DECIMATOR(M, RESPONSE) initializes the filter designer
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
%   D = FDESIGN.DECIMATOR(M, RESPONSE, SPEC) initializes the filter
%   designer 'Specification' property with SPEC.  Use SET(D, 'SPECIFICATION')
%   to get all available specifications for the response RESPONSE. 
%
%   Different design and specification types will have different design
%   methods available. Use DESIGNMETHODS(D) to get a list of design
%   methods available for a given SPEC.
%
%   D = FDESIGN.DECIMATOR(M, RESPONSE, SPEC, SPEC1, SPEC2, ...)
%   initializes the filter designer specifications with SPEC1, SPEC2, etc.
%   Use GET(D, 'DESCRIPTION') for a description of SPEC1, SPEC2, etc.
%
%   By default, all frequency specifications are assumed to be in
%   normalized frequency units. Moreover, all magnitude specifications are
%   assumed to be in dB.
%
%   D = FDESIGN.DECIMATOR(...,Fs) provides the sampling frequency of the
%   signal to be filtered. Fs must be specified as a scalar trailing the
%   other numerical values provided. For this case, Fs is assumed to be in
%   Hz as are all other frequency values provided. Note that you don't
%   change the specification string in this case.
%
%   D = FDESIGN.DECIMATOR(...,MAGUNITS) specifies the units for any
%   magnitude specification given in the constructor. MAGUNITS can be one
%   of the following: 'linear', 'dB', or 'squared'. If this argument is
%   omitted, 'dB' is assumed. Note that the magnitude specifications are
%   always converted and stored in dB regardless of how they were
%   specified.
%
%   Some responses have additional constructor inputs. The syntaxes for
%   these special cases are listed below:
%
%   D = FDESIGN.DECIMATOR(M, 'CIC', DD, SPEC, ...) creates a CIC decimator
%   filter designer with the 'DecimationFactor' property set to M, and the
%   'DifferentialDelay' property set to DD.
%
%   D = FDESIGN.DECIMATOR(M, 'CIC compensator', DD, N, R, SPEC, ...)
%   creates a CIC compensator decimator filter designer with the
%   'DecimationFactor' property set to M, the 'DifferentialDelay' property
%   set to DD, the 'NumberOfSections' property set to N, and the
%   'CICRateChangeFactor' property set to R. In general, the design method
%   implements a filter with a passband response equal to an inverse
%   Dirichlet sinc that matches exactly the inverse passband response of a
%   CIC filter with a differential delay equal to DD, number of sections
%   equal to N, and rate change factor equal to R. When R is not specified
%   or is set to 1, the design method implements a filter with a passband
%   response equal to an inverse sinc that is an approximation to the true
%   inverse passband response of the CIC filter.
%
%   D = FDESIGN.DECIMATOR(M, PULSESHAPERESP, SPS, SPEC, ...) where
%   PULSESHAPERESP equals 'Raised Cosine', 'Square Root Raised Cosine', or
%   'Gaussian', creates a pulse shaping decimator filter designer with the
%   'DecimationFactor' property set to M, and the 'SamplesPerSymbol'
%   property set to SPS.
%
%   % Example #1 - Design a CIC decimator for a signal sampled at 19200 Hz
%   % with a differential delay of one and that attenuates alias below 50 Hz
%   % by at least 80 dB.
%   DD  = 1;     % Differential delay
%   Fp  = 50;    % Passband of interest
%   Ast = 80;    % Minimum attenuation of alias components in passband
%   Fs  = 19200; % Sampling frequency for input signal
%   M   = 64;    % Decimation factor
%   d   = fdesign.decimator(M,'cic',DD,'Fp,Ast',Fp,Ast,Fs);
%   hm  = design(d);
%
%   % Example #2 - Design a minimum-order CIC compensator that decimates by
%   % 4 and compensates for the droop in the passband for the CIC from the
%   % previous example. The stopband of the compensator decays like (1/f)^2.
%   Nsecs = hm.NumberOfSections;
%   R     = hm.DecimationFactor;
%   d     = fdesign.decimator(4,'CIC compensator',DD,Nsecs,R,'Fp,Fst,Ap,Ast',50,100,0.1,80,Fs/R);
%   hmc   = design(d, 'equiripple', 'StopbandShape', '1/f', 'StopbandDecay', 2);
%   hmc.Arithmetic = 'fixed';
%   % Analyze filters individually plus compound response
%   hfvt = fvtool(hm,hmc,cascade(hm,hmc),'Fs',[Fs,Fs/M,Fs],'ShowReference','off'); 
%   legend(hfvt,'CIC decimator','CIC compensator','Overall response');
%
%   % Example #3 - Design a minimum-order Nyquist decimator using a Kaiser
%   % window. Compare to a multistage design. 
%   M   = 15;   % Decimation factor. Also the Nyquist band.
%   TW  = 0.05; % Normalized transition width
%   Ast = 40;   % Minimum stopband attenuation in dB
%   d = fdesign.decimator(M,'nyquist',M,TW,Ast);
%   hm  = design(d,'kaiserwin');  
%   hm2 = design(d,'multistage');
%   hfvt = fvtool(hm,hm2); legend(hfvt,'Kaiser window','Multistage')
%
%   % Example #4 - Design a lowpass decimator for a decimation factor of 8.
%   % Compare a single-stage equiripple design with multistage designs.
%   M = 8; % Decimation factor;
%   d = fdesign.decimator(M,'lowpass');
%   hm(1) = design(d,'equiripple');
%   hm(2) = design(d,'multistage','Usehalfbands',true); % Use halfband filters whenever possible
%   hm(3) = design(d,'multistage','Usehalfbands',true,...
%    'HalfbandDesignMethod','iirlinphase'); % Use quasi linear-phase IIR halfbands
%   hfvt = fvtool(hm);
%   legend(hfvt,'Single-stage equiripple','Multistage','Multistage with IIR halfbands')
%
%   See also FDESIGN, FDESIGN/SETSPECS, FDESIGN/DESIGN,
%   FDESIGN/INTERPOLATOR, FDESIGN/RSRC.

%   Copyright 2005-2011 The MathWorks, Inc.

this = fdesign.decimator;

set(this, 'MultirateType', 'Decimator');

needsDefaults = true;

if nargin > 0
    set(this, 'DecimationFactor', M);
    if nargin > 1
        if isa(DT, 'fdesign.pulseshaping')
            error(message('signal:fdesign:decimator:decimator:unsupportedResponse', class( DT )));
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
    multiratedefaults(this.CurrentFDesign, this.DecimationFactor);
end

setspecs(this, varargin{:});

% [EOF]
