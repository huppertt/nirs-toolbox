function C = wcov (X, w, biased)
    % WCOV Weighted covariance matrix. No NaNs allowed.
    % X is a NxD matrix with data points in the rows.
    % w is a vector of length N of nonnegative reals.
    % biased is a logical scalar, default=false;
    %
    % WCOV(X,W) or WCOV(X,W,0) returns the unbiased weighted covariance
    % matrix.
    % WCOV(X,W,1) returns the biased (maximum likelihood) weighted covariance
    % matrix.
    
    %  Copyright 2014 The MathWorks, Inc.

    if nargin<3
        biased = false;               	% Default is unbiased.
    end
    w = w(:)/sum(w);                  	% make w a normalized column vec.
    Y = bsxfun(@minus, X, w'*X);        % deviations from weighted mean.
    C = Y'*bsxfun(@times, Y, w);        % Weighted Cov.
    % Maybe correct for bias
    if ~biased
        C = C/(1-w'*w);                 % w'*w is sum of squared weights.
    end
    
    % Set entries for constant columns exactly to 0
    constCols = ~var(X);
    C(constCols,:) = 0;
    C(:,constCols) = 0;
    
    % Make it exactly symmetric. Necessary because the weights cause
    % numerical errors.
    C = tril(C) + tril(C,-1)';
end


% REFERENCES

% An approximation to the unbiased weighted sample covariance was used: 
% http://en.wikipedia.org/wiki/Weighted_mean#Weighted_sample_variance.
% This is the same equation used in classreg.learning.internal.wnanvar,
% and the original reference given is:

% Mark Galassi, Jim Davies, James Theiler, Brian Gough, Gerard Jungman,
% Michael Booth, and Fabrice Rossi. GNU Scientific Library - Reference
% manual, Version 1.15, 2011. Sec. 21.7 Weighted Samples.
% http://www.gnu.org/software/gsl/manual/html_node/Weighted-Samples.html

