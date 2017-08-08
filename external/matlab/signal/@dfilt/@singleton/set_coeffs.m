function c = set_coeffs(this,c)                                               
%SET_COEFFS Set the coefficients.                                          

%   Author(s): V. Pellissier                                          
%   Copyright 1988-2005 The MathWorks, Inc.


error(nargchk(2,2,nargin,'struct'));                                                

% Always store as a row                                                    
c = c(:).';                                                                

clearmetadata(this);

% [EOF]
