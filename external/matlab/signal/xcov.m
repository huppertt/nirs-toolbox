function [xycov,lags] = xcov(x,y,option1,option2)
%XCOV   Cross-covariance function estimates.
%   XCOV(A,B), where A and B are length M vectors, returns the
%   length 2*M-1 cross-covariance sequence in a column vector.
%
%   XCOV(A), when A is a vector, is the auto-covariance sequence.
%   XCOV(A), when A is an M-by-N matrix, is a large matrix with
%   2*M-1 rows whose N^2 columns contain the cross-covariance
%   sequences for all combinations of the columns of A.
%   The zeroth lag of the output covariance is in the middle of the 
%   sequence, at element or row M.
%
%   The cross-covariance is the cross-correlation function of
%   two sequences with their means removed:
%        C(m) = E[(A(n+m)-MA)*conj(B(n)-MB)] 
%   where MA and MB are the means of A and B respectively.
%
%   XCOV(...,MAXLAG) computes the (auto/cross) covariance over the
%   range of lags:  -MAXLAG to MAXLAG, i.e., 2*MAXLAG+1 lags.
%   If missing, default is MAXLAG = M-1.
%
%   [C,LAGS] = XCOV(...) returns a vector of lag indices (LAGS).
%   
%   XCOV(...,SCALEOPT), normalizes the covariance according to SCALEOPT:
%       biased   - scales the raw cross-covariance by 1/M.
%       unbiased - scales the raw covariance by 1/(M-abs(k)), where k 
%                  is the index into the result.
%       coeff    - normalizes the sequence so that the covariances at 
%                  zero lag are identically 1.0.
%       none     - no scaling (this is the default).
%
%   % Example:
%   %   Display the estimated autocovariance of white Gaussian noise cww(m)
%   %   for -10 <= m <= 10.
%
%   ww = randn(1000,1);                     % White Gaussian noise
%   [cov_ww,lags] = xcov(ww,10,'coeff');    % Cross-covariance
%   stem(lags,cov_ww)                   
%   
%   See also XCORR, CORRCOEF, CONV, COV and XCORR2.

%   Author(s): L. Shure, 1-9-88
%   Copyright 1988-2012 The MathWorks, Inc.

%   References:
%     [1] J.S. Bendat and A.G. Piersol, "Random Data:
%         Analysis and Measurement Procedures", John Wiley
%         and Sons, 1971, p.332.
%     [2] A.V. Oppenheim and R.W. Schafer, Digital Signal 
%         Processing, Prentice-Hall, 1975, pg 539.


if nargin == 1
    [xycov,l] = xcorr( bsxfun(@minus, x, mean(x)) );
elseif nargin == 2
    if ischar(y)||(~ischar(y)&&length(y)==1)
        [xycov,l] = xcorr(bsxfun(@minus, x, mean(x)), y);
    else
        [xycov,l] = xcorr(bsxfun(@minus, x, mean(x)), ...
                bsxfun(@minus, y, mean(y)) );
    end
elseif nargin == 3
    if length(y)==1
        [xycov,l] = xcorr(bsxfun(@minus, x, mean(x)), y, option1);
    else
        [xycov,l] = xcorr(bsxfun(@minus, x, mean(x)), ...
                bsxfun(@minus, y, mean(y)), option1);
    end
elseif nargin == 4
    [xycov,l] = xcorr(bsxfun(@minus, x, mean(x)), ...
                bsxfun(@minus, y, mean(y)), option1, option2);
end
if nargout > 1
    lags = l;
end

