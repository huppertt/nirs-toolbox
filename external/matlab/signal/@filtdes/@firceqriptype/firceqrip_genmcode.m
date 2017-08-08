function [params, values, descs, args] = firceqrip_genmcode(h, d)
%FIRCEQRIP_GENMCODE

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

% Get frequency spec, it has been prenormalized
freqSpecType = get(d,'freqSpecType');
propname     = determine_dynamicprop(d,freqSpecType,set(d,'freqSpecType'));

params = {'N', propname, 'slope'};
values = {getmcode(d, 'order'), getmcode(d, propname), getmcode(d, 'stopbandslope')};
descs  = {'', '', 'Stopband Slope'};

args = sprintf('N, %s%s, %%s, ''slope'', slope', propname, getfsstr(d));

if ~strcmpi(freqSpecType,'cutoff'),
    args = sprintf('%s, spectype', args);
    params = {params{:}, 'spectype'};
    values = {values{:}, sprintf('''%s''', freqSpecType)};
    descs  = {descs{:}, 'Frequency Specification Type'};
end

if strcmpi(get(d,'minPhase'),'on'),
    args = sprintf('%s, ''min''', args);
end

% [EOF]
