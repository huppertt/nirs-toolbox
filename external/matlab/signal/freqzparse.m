function [NorW, UC, Fs, FsS] = freqzparse(varargin)
%FREQZPARSE parse the inputs to freqz
%   [NorW, UC, Fs, FsS] = freqzparse(varargin{:}) Returns Nfft or the W
%   vector in the first output, UC ('whole' vs 'half') in the 2nd, the Fs
%   in the 3rd and a boolean flag that is true when Fs is specified.  The
%   default values are 8192, 'half' and [] respectively.

%   Copyright 1988-2013 The MathWorks, Inc.

NorW = 8192;
UC   = 'half';
Fs   = [];

% FsS is for freqz.  when you specify [] in freqz it defaults to an Fs of 1
% instead of normalizing.  so if the user did not specify an fs we need to
% pass that information back to the caller of this function.
FsS  = false;

switch nargin
    case 1,
        
        % If there is only 1 input it either NorW or the unitcircle
        % specification.
        if isnumeric(varargin{1}),
            NorW = varargin{1};
        elseif ischar(varargin{1}),
            UC = varargin{1};
        end
    case 2,
        
        % If there are 2 inputs the first is NorW and the 2nd is either the
        % Fs or the unit circle specification
        NorW = varargin{1};
        if isnumeric(varargin{2}),
            Fs  = varargin{2};
            FsS = true;
        elseif ischar(varargin{2}),
            UC = varargin{2};
        end
    case 3,
        
        % If there are 3 inputs they are exactly known (NorW, UC, Fs)
        NorW = varargin{1};
        FsS  = true;
        if ischar(varargin{2}),
            UC = varargin{2};
            Fs = varargin{3};
        elseif ischar(varargin{3}),
            UC = varargin{3};
            Fs = varargin{2};
        else
            error(message('signal:freqzparse:InvalidFrequencyRange'));
        end
end

UC = validatestring(UC,{'half','whole'});

% [EOF]
