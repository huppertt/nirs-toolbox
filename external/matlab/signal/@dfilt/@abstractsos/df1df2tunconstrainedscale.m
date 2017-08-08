function df1df2tunconstrainedscale(Hd,opts,L)
%DF1DF2TUNCONSTRAINEDSCALE   

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

% Compute unconstrained scaling

% Calculate cumulative scaling factors
c = cscalefactors(Hd,opts);
c = c(:);
% Prepend a cumulative scale factor of 1
c = [1; c];

% Compute unconstrained scaling using scale factors
scalefact = c(2:L+1)./c(1:L);
if strcmpi(opts.ScaleValueConstraint,'unit'),
    scalefact = repmat(scalefact,1,3);
    Hd.sosMatrix(1:L,1:3) = Hd.sosMatrix(1:L,1:3).*scalefact;
    Hd.ScaleValues(L+1) = Hd.ScaleValues(L+1)/c(L+1);

    % Incorporate existing scale values
    scalevals = repmat(Hd.ScaleValues(1:L),1,3);
    Hd.sosMatrix(1:L,1:3) = Hd.sosMatrix(1:L,1:3).*scalevals;
    Hd.ScaleValues(1:L) = ones(L,1);
else
    Hd.ScaleValues(1:L) = Hd.ScaleValues(1:L).*scalefact;
    Hd.ScaleValues(L+1) = Hd.ScaleValues(L+1)/c(L+1);
end


% [EOF]
