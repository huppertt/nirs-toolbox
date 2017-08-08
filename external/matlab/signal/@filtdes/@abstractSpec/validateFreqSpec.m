function varargout = validateFreqSpec(~, d, varargin)
%VALIDATEFREQSPEC Validate the freqspec

%   Copyright 2012 The MathWorks, Inc.

success   = true;
exception = MException.empty;

% Calculate the nyquest frequency
if d.isnormalized
    nyq = 1;
else
    nyq = d.Fs/2;
end

% Go through all the properties and make sure they are all below the
% nyquist frequency.
for indx = 1:numel(varargin)
    try
        value = d.(varargin{indx});
    catch ME
        switch ME.identifier
            case {'MATLAB:noSuchMethodOrField', ...
                'MATLAB:class:GetDenied',...
                'MATLAB:UndefinedFunction'}
                continue;
            otherwise
                rethrow(ME);
        end
    end
    if value > nyq
        success = false;
        exception = MException('signal:fdatool:InvalidFreqSpec', ...
            fdatoolmessage('InvalidFreqSpec', varargin{indx}, num2str(nyq)));
    end
end

if nargout
    varargout = {success, exception};
elseif ~success
    throw(exception);
end

% [EOF]
