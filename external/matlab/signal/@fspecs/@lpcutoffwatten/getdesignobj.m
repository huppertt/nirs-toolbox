function designobj = getdesignobj(~, str, sigonlyflag) 
%GETDESIGNOBJ   Get the design object.

%   Copyright 1988-2013 The MathWorks, Inc.

if nargin < 2
    str = [];
    sigonlyflag = false;
elseif nargin < 3
    sigonlyflag = false;
end

if isfdtbxinstalled && ~sigonlyflag
    %#function fdfmethod.eqriplpcutoff
    designobj.equiripple = 'fdfmethod.eqriplpcutoff';
else
    designobj = [];
end
%#function fmethod.firclslp
designobj.fircls = 'fmethod.firclslp';

if ~isempty(str)
    designobj = designobj.(str);
end

% [EOF]
