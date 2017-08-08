function zi = filtfilt_initial_conditions(b,a)
%FILTFILT_INITIAL_CONDITIONS(B,A) returns the initial conditions for FILTFILT.

%   Copyright 1988-2009 The MathWorks, Inc.

% set up filter's initial conditions to remove dc offset problems at the 
% beginning and end of the sequence

nb = length(b);
na = length(a);
nfilt = max(nb,na);
% convert to row and zero-pad if necessary
    b = [b(:).', zeros(1,nfilt-nb)];
    a = [a(:).', zeros(1,nfilt-na)];
    
% use sparse matrix to solve system of linear equations for initial conditions
% zi are the steady-state states of the filter b(z)/a(z) in the state-space 
% implementation of the 'filter' command.
    rows = [1:nfilt-1  2:nfilt-1  1:nfilt-2];
    cols = [ones(1,nfilt-1) 2:nfilt-1  2:nfilt-1];
    data = [1+a(2) a(3:nfilt) ones(1,nfilt-2)  -ones(1,nfilt-2)];
    sp = sparse(rows,cols,data);
    zi = full(sp \ ( b(2:nfilt).' - a(2:nfilt).'*b(1) ));
% non-sparse:
%  zi = ( eye(nfilt-1) - [-a(2:nfilt).' [eye(nfilt-2); zeros(1,nfilt-2)]] ) \ ...
%       ( b(2:nfilt).' - a(2:nfilt).'*b(1) );
