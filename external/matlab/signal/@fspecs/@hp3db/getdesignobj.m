function designobj = getdesignobj(~, str, sigonlyflag)
%GETDESIGNOBJ   Get the design object.

%   Copyright 2005-2013 The MathWorks, Inc.

if nargin < 2
    str = [];
    sigonlyflag = false;
elseif nargin < 3
    sigonlyflag = false;
end

%#function fmethod.butterhp
designobj.butter = 'fmethod.butterhp';

if isfdtbxinstalled && ~sigonlyflag
    %#function fdfmethod.firmaxflathp
    designobj.maxflat = 'fdfmethod.firmaxflathp';
end

if ~isempty(str)
    designobj = designobj.(str);
end

% [EOF]
