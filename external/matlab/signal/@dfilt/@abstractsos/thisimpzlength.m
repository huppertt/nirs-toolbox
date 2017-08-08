function len = thisimpzlength(Hd, varargin)
%THISIMPZLENGTH Length of the impulse response for a digital filter.
%   THISIMPZLENGTH(Hd) returns the length of the impulse response of 
%   the filter defined by Hd.
%  
%   THISIMPZLENGTH(Hd,TOL) will specify the tolerance for greater or 
%   less accuracy.  By default, TOL = 5e-5.
%  
%   See also IMPZ.

%   Author: R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

error(nargchk(1,2,nargin,'struct'));

warnsv(Hd);

% Initialize length
firlen=1;
iirlen=1;

% Convert the filter to a transfer function.
for k=1:nsections(Hd)
    
    % Get the transfer function coefficients
    b=Hd.sosMatrix(k,1:3);
    a=Hd.sosMatrix(k,4:6);
    
    if signalpolyutils('isfir',b,a),
        % Add the length of each FIR section
        firlen = firlen + length(b) - 1;
    else 
        
        % Keep the maximum length of all IIR sections
        iirlen = max(iirlen, impzlength(b,a,varargin{:}));
    end
end

% Use the longest of FIR or IIR
len=max(firlen,iirlen);
