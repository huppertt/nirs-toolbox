function p = propstosync(this)
%PROPSTOSYNC   

%   Copyright 2005-2012 The MathWorks, Inc.

p = get(classhandle(this), 'Properties');

p = find(p, 'AccessFlags.PublicSet', 'On', '-not', 'Name', 'Fs',...
  '-not','Name','NormalizedFrequency','-not','Name','FromFilterDesigner');

p = get(p, 'Name');

if ~iscellstr(p)
    p = {p};
end

p = thispropstosync(this,p);

