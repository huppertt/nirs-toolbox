function len = thisimpzlength(Hd, varargin)
%THISIMPZLENGTH Length of the impulse response for a digital filter.
%   IMPZLENGTH(Hd) returns the length of the impulse response of 
%   the filter defined by Hd.
%  
%   IMPZLENGTH(Hd,TOL) will specify the tolerance for greater or 
%   less accuracy.  By default, TOL = 5e-5.
%  
%   See also IMPZ.

%   Author: R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.

% This should be private

error(nargchk(1,2,nargin,'struct'));

% Initialize length
firlen = 1;
iirlen = 1;


for k=1:length(Hd.Stage)
    
    % Convert the filter to a transfer function.
    [b, a] = tf(Hd.Stage(k));
    
    if signalpolyutils('isfir',b,a),
        % Add the length of each FIR section
        firlen = firlen + length(b) - 1;
    else
        
        % Keep the maximum length of all IIR sections
        iirlen = max(iirlen,impzlength(b,a,varargin{:}));
    end
end

    
% Use the longest of FIR or IIR
len = max(firlen,iirlen);



