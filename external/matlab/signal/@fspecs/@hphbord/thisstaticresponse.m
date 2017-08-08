function thisstaticresponse(this)
%THISSTATICRESPONSE <short description>

%   Copyright 2007 The MathWorks, Inc.

if this.NormalizedFrequency, str = '.5';
else,                        str = 'Fs/4'; end

staticrespengine('drawpassband',   hax, [.55 1], [.9 1.1]);
staticrespengine('drawstopband',   hax, [0 .45]);
staticrespengine('drawfreqlabels', hax, .5, str);


% [EOF]
