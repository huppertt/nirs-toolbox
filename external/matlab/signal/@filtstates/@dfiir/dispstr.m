function s = dispstr(Hs,offset)
%DISPSTR Display string of states.
%   DISPSTR(Hs) returns a string that can be used to display the states
%   of a Filter States object (FILTSTATES.DFIIR).
  
%   Author: P. Costa
%   Copyright 1988-2004 The MathWorks, Inc.
  
[NumStr,DenStr] = getstr(Hs);
s = ['Numerator:',blanks(2),NumStr,'\n',blanks(offset),'Denominator:',DenStr];

% -------------------------------------------------------------------------
function [NumStr,DenStr] = getstr(h)

% Left two separate checks since the FILTSTATES object (when not contained
% in a DFILT) can have one property as a FI and the other as a non-FI.
if strcmpi(class(h.Numerator),'embedded.fi'), 
    ncls = 'fi'; 
else
    ncls = class(h.Numerator);
end

if strcmpi(class(h.Denominator),'embedded.fi'), 
    dcls = 'fi'; 
else
    dcls = class(h.Denominator);
end


szNumS = size(h.Numerator);
szDenS = size(h.Denominator);
NumStr = ['[',num2str(szNumS(1)),'x',num2str(szNumS(2)),' ',ncls,']'];
DenStr = ['[',num2str(szDenS(1)),'x',num2str(szDenS(2)),' ',dcls,']'];
         