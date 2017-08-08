function fr = dmom_whichframes(h)
%WHICHFRAMES  Return constructors of frames needed for FDATool.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Call super's method
fr = super_whichframes(h);

% Change isMinOrd in filter order frame
indx = find(strcmpi({fr.constructor},'siggui.filterorder'));
if isdynpropenab(h,'orderMode');
    fr(indx).setops  = {'isMinOrd',1};
end
