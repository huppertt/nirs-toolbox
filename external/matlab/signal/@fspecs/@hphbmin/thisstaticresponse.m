function thisstaticresponse(this, hax)
%THISSTATICRESPONSE   

%   Copyright 2007 The MathWorks, Inc.

magunits = getappdata(hax, 'magunits');
if strcmpi(magunits, 'linear'), inputs = {'D_{stop}', 'right'};
else,                           inputs = {'A_{stop}'}; end

if this.NormalizedFrequency, str = '.5';
else,                        str = 'Fs/4'; end

staticrespengine('drawstopband',   hax, [0 .45], inputs{:});
staticrespengine('drawtransition', hax, [.45 .55]);
staticrespengine('drawpassband',   hax, [.55 1], [.9 1.1]);
staticrespengine('drawfreqlabels', hax, .5, str);


% [EOF]
