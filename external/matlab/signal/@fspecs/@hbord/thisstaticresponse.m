function thisstaticresponse(this, hax)
%THISSTATICRESPONSE   

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

if this.NormalizedFrequency, str = '.5';
else,                        str = 'Fs/4'; end

staticrespengine('drawpassband',   hax, [0   .45], [.9 1.1]);
staticrespengine('drawstopband',   hax, [.55 1]);
staticrespengine('drawfreqlabels', hax, .5, str);

% [EOF]
