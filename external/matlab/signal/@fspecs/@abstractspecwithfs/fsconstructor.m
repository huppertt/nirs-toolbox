function fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin)
%FSCONSTRUCTOR   Base constructor for all specs with Fs.

%   Author(s): R. Losada
%   Copyright 2003-2011 The MathWorks, Inc.
%     


error(nargchk(0,nargsnoFs+1,length(varargin),'struct'));


freqargs = [varargin{fstart:min(length(varargin),fstop)}];

if length(varargin) > nargsnoFs
    % Fs specified
    Fs = varargin{nargsnoFs+1};
    if any(freqargs > Fs/2),
        error(message('signal:fspecs:abstractspecwithfs:fsconstructor:SpecBeyondNyquist'));
    end
else
    % Fs not specified
    if ~isempty(freqargs) && any(freqargs > 1),
        error(message('signal:fspecs:abstractspecwithfs:fsconstructor:invalidSpec'));
    end
end

this.ResponseType = respstr;

this.setspecs(varargin{:});

% [EOF]
