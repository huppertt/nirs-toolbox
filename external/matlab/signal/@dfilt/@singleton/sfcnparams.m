function params = sfcnparams(Hd, library)
%SFCNPARAMS Returns the parameters for SDSPFILTER

% Author(s): J. Schickler
% Copyright 1988-2005 The MathWorks, Inc.

%% sdspfilter requires the following 8 input arguments:
%% (Please refer to sdspfilter.c for more information)
%% CoeffsFromMask:   1 or 0.  Set this to 1 because the DFD block
%%                   does not have time-varying coeffs.
%% FilterType:       Set this to 0 to indicate an IIR filter
%% FilterStruct:     Set this to 5 to indicate a DF2T structure
%% FilterUpdateRate: Set this to 0 to indicate one filt per frame
%% CoeffsNorm:       Set this to 1 if coefficients are normalized
%% Coeff1:           variable name for the first coeff set
%% Coeff2:           variable name for the second coeff set
%% IC:               0 because this isn't supported yet

if nargin < 2,
    library = 'dsparch4';
end

% Get the filter structure specific information from the subclasses.
[filtertype, filterstrt, num, den, ic] = thissfcnparams(Hd);

lib = 1;
if ~isempty(strfind(library, 'dsparch3')),
    lib = 0;
end

params = sprintf('1, %d, %d, 0, 1, [%s], [%s], %s, %d', ...
    filtertype, filterstrt, num, den, ic, lib);

% [EOF]
