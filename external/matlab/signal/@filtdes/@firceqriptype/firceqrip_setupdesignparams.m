function args = firceqrip_setupdesignparams(h,d)
%FIRCEQRIP_SETUPDESIGNPARAMS

%   Author(s): R. Losada, J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

% Set up design params
N = get(d,'order');

% Get frequency spec, it has been prenormalized
freqSpecType = get(d,'freqSpecType');
propname = determine_dynamicprop(d,freqSpecType,set(d,'freqSpecType'));
Fspec = get(d,propname);

slope = get(d,'stopbandSlope');

args = {N,Fspec,'replaceMe','slope',slope};

if ~strcmpi(freqSpecType,'cutoff'),
    args = {args{:},freqSpecType};
end
    
if strcmpi(get(d,'minPhase'),'on'),
    args = {args{:},'min'};
end
