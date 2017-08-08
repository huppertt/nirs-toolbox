function strs = whichspecobjs(h)
%WHICHSPECOBJS Determine which specs objects are used by this class.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

strs = firceqrip_whichspecobjs(h);
strs = {strs{:},'filtdes.hpmagfir'};
