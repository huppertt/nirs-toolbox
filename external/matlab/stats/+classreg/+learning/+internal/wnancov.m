function C = wnancov (X, w, biased)
    % WNANCOV Weighted covariance matrix.
    % X is a NxD matrix with data points in the rows.
    % w is a vector of length N of nonnegative reals.
    % biased is a logical scalar, default=false;
    %
    % C = WNANCOV(X,W) or C = WNANCOV(X,W,0) returns the unbiased weighted
    % covariance matrix (Default).
    %
    % C = WNANCOV(X,W,1) returns the biased (maximum likelihood) weighted
    % covariance matrix.

    %  Copyright 2014 The MathWorks, Inc.

    if nargin<3
        biased = false;                     % Default is unbiased.
    end
    % Delete rows containing NaNs in X or w, then call wcov.
    nanrows = any(isnan(X),2) | isnan(w);
    if any(nanrows)
        X = X(~nanrows,:);
        w = w(~nanrows);
    end
    C = classreg.learning.internal.wcov(X, w, biased);
end

