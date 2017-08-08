function designobj = getdesignobj(~, str, sigonlyflag)
%GETDESIGNOBJ   Get the designobj.

%   Copyright 1988-2013 The MathWorks, Inc.

if nargin < 2
    str = [];
    sigonlyflag = false;
elseif nargin < 3
    sigonlyflag = false;
end
    
%#function fmethod.butterlpmin
%#function fmethod.cheby1lpmin
%#function fmethod.cheby2lpmin
%#function fmethod.elliplpmin
designobj.butter     = 'fmethod.butterlpmin';
designobj.cheby1     = 'fmethod.cheby1lpmin';
designobj.cheby2     = 'fmethod.cheby2lpmin';
designobj.ellip      = 'fmethod.elliplpmin';

if isfdtbxinstalled && ~sigonlyflag
    %#function fdfmethod.eqriplpmin
    %#function fdfmethod.ifirlpmin
    %#function fdfmethod.multistage
    designobj.equiripple = 'fdfmethod.eqriplpmin';
    designobj.ifir       = 'fdfmethod.ifirlpmin';
    designobj.multistage = 'fdfmethod.multistage';
else
    %#function fmethod.eqriplpmin
    designobj.equiripple = 'fmethod.eqriplpmin';
end
%#function fmethod.kaiserlpmin
designobj.kaiserwin  = 'fmethod.kaiserlpmin';

if ~isempty(str)
    designobj = designobj.(str);
end

% [EOF]
