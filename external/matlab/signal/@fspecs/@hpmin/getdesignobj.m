function designobj = getdesignobj(~, str, sigonlyflag)
%GETDESIGNOBJ   Get the designobj.

%   Copyright 1988-2013 The MathWorks, Inc.

if nargin < 2
    str = [];
    sigonlyflag = false;
elseif nargin < 3
    sigonlyflag = false;
end

%#function fmethod.butterhpmin
%#function fmethod.cheby1hpmin
%#function fmethod.cheby2hpmin
%#function fmethod.elliphpmin
designobj.butter     = 'fmethod.butterhpmin';
designobj.cheby1     = 'fmethod.cheby1hpmin';
designobj.cheby2     = 'fmethod.cheby2hpmin';
designobj.ellip      = 'fmethod.elliphpmin';

if isfdtbxinstalled && ~sigonlyflag
    %#function fdfmethod.eqriphpmin
    %#function fdfmethod.ifirhpmin
    designobj.equiripple = 'fdfmethod.eqriphpmin';
    designobj.ifir       = 'fdfmethod.ifirhpmin';
else
    %#function fmethod.eqriphpmin
    designobj.equiripple = 'fmethod.eqriphpmin';
end
%#function fmethod.kaiserhpmin
designobj.kaiserwin  = 'fmethod.kaiserhpmin';

if ~isempty(str)
    designobj = designobj.(str);
end

% [EOF]
