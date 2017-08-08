function designobj = getdesignobj(this, str)
%GETDESIGNOBJ   Get the designobj.

%   Copyright 2007 The MathWorks, Inc.

if isfdtbxinstalled
    
    %#function fdfmethod.windowhphbord
    %#function fdfmethod.butterhphalfband
    designobj.window = 'fdfmethod.windowhphbord';
    designobj.butter = 'fdfmethod.butterhphalfband';
else
    designobj = [];
end

if nargin > 1
    designobj = designobj.(str);
end

% [EOF]