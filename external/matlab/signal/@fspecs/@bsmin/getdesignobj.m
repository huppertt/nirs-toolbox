function designobj = getdesignobj(~, str, sigonlyflag)
%GETDESIGNOBJ   Get the designobj.

%   Copyright 1988-2013 The MathWorks, Inc.

if nargin < 2
    str = [];
    sigonlyflag = false;
elseif nargin < 3
    sigonlyflag = false;
end

%#function fmethod.butterbsmin
%#function fmethod.cheby1bsmin
%#function fmethod.cheby2bsmin
%#function fmethod.ellipbsmin
designobj.butter     = 'fmethod.butterbsmin';
designobj.cheby1     = 'fmethod.cheby1bsmin';
designobj.cheby2     = 'fmethod.cheby2bsmin';
designobj.ellip      = 'fmethod.ellipbsmin';

if isfdtbxinstalled && ~sigonlyflag
    %#function fdfmethod.eqripbsmin
    designobj.equiripple = 'fdfmethod.eqripbsmin';
else
    %#function fmethod.eqripbsmin
    designobj.equiripple = 'fmethod.eqripbsmin';
end


%#function fmethod.kaiserbsmin
designobj.kaiserwin  = 'fmethod.kaiserbsmin';

if ~isempty(str)
    designobj = designobj.(str);
end

% [EOF]
