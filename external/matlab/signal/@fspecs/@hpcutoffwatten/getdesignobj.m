function designobj = getdesignobj(~, str, sigonlyflag)
%GETDESIGNOBJ   Get the design object.

%   Copyright 1988-2013 The MathWorks, Inc.

if nargin < 2
    str = [];
    sigonlyflag = false;
elseif nargin < 3
    sigonlyflag = false;
end

if isfdtbxinstalled  && ~sigonlyflag   
    %#function fdfmethod.eqriphpcutoff
    designobj.equiripple = 'fdfmethod.eqriphpcutoff';
else
    designobj = [];
end
%#function fmethod.firclshp
designobj.fircls = 'fmethod.firclshp';

if ~isempty(str)
    designobj = designobj.(str);
end

% [EOF]
