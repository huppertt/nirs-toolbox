function [isvalid,errmsg,msgid] = validate(h)
%VALIDATE   Validate specs.

%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.

% Populate defaults

if h.NormalizedFrequency,
    upperlimit = 1;
else
    % Fs is used 
    upperlimit = h.Fs/2;
end

% Get all frequency specs
fspecs = get(h,props2normalize(h));

% Make it a vector
fspecs = [fspecs{:}];

if any(fspecs > upperlimit),
    error(message('signal:fspecs:abstractspecwithfs:validate:invalidSpec', sprintf( '%0.5g', upperlimit )));
end

[isvalid, errmsg, msgid] = thisvalidate(h);

% [EOF]
