function thisstaticresponse(this, hax)
%THISSTATICRESPONSE   

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

magunits = getappdata(hax, 'magunits');
if strcmpi(magunits, 'linear'), inputs = {'D_{stop}', 'right'};
else,                           inputs = {'A_{stop}'}; end

if this.NormalizedFrequency, str = '.5';
else,                        str = 'Fs/4'; end

staticrespengine('drawpassband',   hax, [0   .45], [.9 1.1]);
staticrespengine('drawtransition', hax, [.45 .55]);
staticrespengine('drawstopband',   hax, [.55 1], inputs{:});
staticrespengine('drawfreqlabels', hax, .5, str);

% [EOF]
