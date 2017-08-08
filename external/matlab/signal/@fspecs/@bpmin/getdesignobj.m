function designobj = getdesignobj(~, str, sigonlyflag)
%GETDESIGNOBJ   Get the designobj.

%   Copyright 1988-2013 The MathWorks, Inc.

if nargin < 2
    str = [];
    sigonlyflag = false;
elseif nargin < 3
    sigonlyflag = false;
end

%#function fmethod.butterbpmin
%#function fmethod.cheby1bpmin
%#function fmethod.cheby2bpmin
%#function fmethod.ellipbpmin
designobj.butter     = 'fmethod.butterbpmin';
designobj.cheby1     = 'fmethod.cheby1bpmin';
designobj.cheby2     = 'fmethod.cheby2bpmin';
designobj.ellip      = 'fmethod.ellipbpmin';

if isfdtbxinstalled && ~sigonlyflag
    %#function fdfmethod.eqripbpmin
    designobj.equiripple = 'fdfmethod.eqripbpmin';
else
    %#function fmethod.eqripbpmin
    designobj.equiripple = 'fmethod.eqripbpmin';
end

%#function fmethod.kaiserbpmin
designobj.kaiserwin  = 'fmethod.kaiserbpmin';

if ~isempty(str)
    designobj = designobj.(str);
end

% [EOF]
