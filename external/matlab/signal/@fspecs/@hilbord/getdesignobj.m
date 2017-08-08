function designobj = getdesignobj(~, str, sigonlyflag)
%GETDESIGNOBJ   Get the designobj.

%   Copyright 2005-2013 The MathWorks, Inc.

if nargin < 2
    str = [];
    sigonlyflag = false;
elseif nargin < 3
    sigonlyflag = false;
end

%#function fmethod.eqriphilbord
%#function fmethod.firlshilbord
designobj.equiripple = 'fmethod.eqriphilbord';
designobj.firls      = 'fmethod.firlshilbord';

if isfdtbxinstalled && ~sigonlyflag    
    %#function fdfmethod.eqriphilbord
    %#function fdfmethod.elliphilbertfpass
    %#function fdfmethod.iirlinphasehilbertfpass    
    designobj.equiripple =  'fdfmethod.eqriphilbord';
    designobj.ellip       = 'fdfmethod.elliphilbertfpass';
    designobj.iirlinphase = 'fdfmethod.iirlinphasehilbertfpass';    
end

if ~isempty(str)
    designobj = designobj.(str);
end

% [EOF]
