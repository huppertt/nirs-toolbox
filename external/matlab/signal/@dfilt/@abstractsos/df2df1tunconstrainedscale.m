function df2df1tunconstrainedscale(Hd,opts,L)
%DF2DF1TUNCONSTRAINEDSCALE   

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

% Calculate cumulative scaling factors
c = cscalefactors(Hd,opts);
c1 = c(1,1:L).'; % For no scale values, only use top row

% Append a final cumulative scale factor
c1 = [c1;1];
Hd.ScaleValues(1) = Hd.ScaleValues(1)*c1(1);

% Compute unconstrained scaling using scale factors
scalefact = c1(2:L+1)./c1(1:L);

if strcmpi(opts.ScaleValueConstraint,'unit'),
    scalefact = repmat(scalefact,1,3);
    Hd.sosMatrix(1:L,1:3) = Hd.sosMatrix(1:L,1:3).*scalefact;

    % Incorporate existing scale values except the first
    scalevals = repmat(Hd.ScaleValues(2:L),1,3);
    Hd.sosMatrix(1:L-1,1:3) = Hd.sosMatrix(1:L-1,1:3).*scalevals;
    Hd.sosMatrix(L,1:3) = Hd.sosMatrix(L,1:3).*Hd.ScaleValues(L+1);
    Hd.ScaleValues(2:L+1) = ones(L,1);
else
    Hd.ScaleValues(2:L+1) = Hd.ScaleValues(2:L+1).*scalefact;

    c2 = c(2,1:L).'; % Scale to input of scale values as well
    addscalefact = c2./c1(1:L);
    Hd.ScaleValues(2:L+1) = Hd.ScaleValues(2:L+1)./addscalefact;
    addscalefact = repmat(addscalefact,1,3);
    Hd.sosMatrix(1:L,1:3) = Hd.sosMatrix(1:L,1:3).*addscalefact;
end


% [EOF]
