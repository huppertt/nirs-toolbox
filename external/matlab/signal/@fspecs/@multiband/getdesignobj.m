function designobj = getdesignobj(this, str, sigonlyflag)
%GETDESIGNOBJ Get the design object.

%   Copyright 2005-2013 The MathWorks, Inc.

if nargin < 2
    str = [];
    sigonlyflag = false;
elseif nargin < 3
    sigonlyflag = false;
end

%#function fmethod.firlsmultiband
designobj.firls = 'fmethod.firlsmultiband';

if isfdtbxinstalled && ~sigonlyflag
    %#function fdfmethod.eqripmultiband
    designobj.equiripple = 'fdfmethod.eqripmultiband';
    [F, ~] = getmask(this);
    if all(F>=0),
        %#function fdfmethod.lpnormmultiband1
        designobj.iirlpnorm = 'fdfmethod.lpnormmultiband1';
    end
else    
    %#function fmethod.eqripmultiband
    designobj.equiripple = 'fmethod.eqripmultiband';
end

if ~isempty(str)
    designobj = designobj.(str);
end

% [EOF]
