function designobj = getdesignobj(this, str, sigonlyflag)
%GETDESIGNOBJ Get the design object.

%   Copyright 2005-2013 The MathWorks, Inc.

if nargin < 2
    str = [];
    sigonlyflag = false;
elseif nargin < 3
    sigonlyflag = false;
end

%#function fmethod.freqsamparbmag
%#function fmethod.firlssbarbmag
designobj.freqsamp = 'fmethod.freqsamparbmag';
designobj.firls = 'fmethod.firlssbarbmag';

if isfdtbxinstalled && ~sigonlyflag
    %#function fdfmethod.eqripsbarbmag
    designobj.equiripple = 'fdfmethod.eqripsbarbmag';
    if all(this.Frequencies>=0),
        %#function fdfmethod.lpnormsbarbmag1
        designobj.iirlpnorm = 'fdfmethod.lpnormsbarbmag1';
    end
else
    %#function fmethod.eqripsbarbmag
    designobj.equiripple = 'fmethod.eqripsbarbmag';
end

if ~isempty(str)
    designobj = designobj.(str);
end

% [EOF]
