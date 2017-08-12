function y = prctile(x,p,dim)
%PRCTILE Percentiles of a sample.
%   Y = PRCTILE(X,P) returns percentiles of the values in X.  P is a scalar
%   or a vector of percent values.  When X is a vector, Y is the same size
%   as P, and Y(i) contains the P(i)-th percentile.  When X is a matrix,
%   the i-th row of Y contains the P(i)-th percentiles of each column of X.
%   For N-D arrays, PRCTILE operates along the first non-singleton
%   dimension.
%
%   Y = PRCTILE(X,P,DIM) calculates percentiles along dimension DIM.  The
%   DIM'th dimension of Y has length LENGTH(P).
%
%   Percentiles are specified using percentages, from 0 to 100.  For an N
%   element vector X, PRCTILE computes percentiles as follows:
%      1) The sorted values in X are taken as the 100*(0.5/N), 100*(1.5/N),
%         ..., 100*((N-0.5)/N) percentiles.
%      2) Linear interpolation is used to compute percentiles for percent
%         values between 100*(0.5/N) and 100*((N-0.5)/N)
%      3) The minimum or maximum values in X are assigned to percentiles
%         for percent values outside that range.
%
%   PRCTILE treats NaNs as missing values, and removes them.
%
%   Examples:
%      y = prctile(x,50); % the median of x
%      y = prctile(x,[2.5 25 50 75 97.5]); % a useful summary of x
%
%   See also IQR, MEDIAN, NANMEDIAN, QUANTILE.

%   Copyright 1993-2015 The MathWorks, Inc.


if ~isvector(p) || numel(p) == 0 || any(p < 0 | p > 100) || ~isreal(p)
    error(message('stats:prctile:BadPercents'));
end

% Figure out which dimension prctile will work along.
sz = size(x);
if nargin < 3 
    dim = find(sz ~= 1,1);
    if isempty(dim)
        dim = 1; 
    end
    dimArgGiven = false;
else
    % Permute the array so that the requested dimension is the first dim.
    nDimsX = ndims(x);
    perm = [dim:max(nDimsX,dim) 1:dim-1];
    x = permute(x,perm);
    % Pad with ones if dim > ndims.
    if dim > nDimsX
        sz = [sz ones(1,dim-nDimsX)];
    end
    sz = sz(perm);
    dim = 1;
    dimArgGiven = true;
end

% If X is empty, return all NaNs.
if isempty(x)
    if isequal(x,[]) && ~dimArgGiven
        y = nan(size(p),'like',x);
    else
        szout = sz; szout(dim) = numel(p);
        y = nan(szout,'like',x);
    end

else
    % Drop X's leading singleton dims, and combine its trailing dims.  This
    % leaves a matrix, and we can work along columns.
    nrows = sz(dim);
    ncols = numel(x) ./ nrows;
    x = reshape(x, nrows, ncols);

    x = sort(x,1);
    n = sum(~isnan(x), 1); % Number of non-NaN values in each column
    
    % For columns with no valid data, set n=1 to get nan in the result
    n(n==0) = 1;

    % If the number of non-nans in each column is the same, do all cols at once.
    if all(n == n(1))
        n = n(1);
        if isequal(p,50) % make the median fast
            if rem(n,2) % n is odd
                y = x((n+1)/2,:);
            else        % n is even
                y = (x(n/2,:) + x(n/2+1,:))/2;
            end
        else
            y = interpColsSame(x,p,n);
        end

    else
        % Get percentiles of the non-NaN values in each column.
        y = interpColsDiffer(x,p,n);
    end

    % Reshape Y to conform to X's original shape and size.
    szout = sz; szout(dim) = numel(p);
    y = reshape(y,szout);
end
% undo the DIM permutation
if dimArgGiven
     y = ipermute(y,perm);  
end

% If X is a vector, the shape of Y should follow that of P, unless an
% explicit DIM arg was given.
if ~dimArgGiven && isvector(x)
    y = reshape(y,size(p)); 
end


function y = interpColsSame(x, p, n)
%INTERPCOLSSAME An aternative approach of 1-D linear interpolation which is
%   faster than using INTERP1Q and can deal with invalid data so long as
%   all columns have the same number of valid entries (scalar n).

% Make p a column vector. Note that n is assumed to be scalar.
if isrow(p)
    p = p';
end

% Form the vector of index values (numel(p) x 1)
r = (p/100)*n;
k = floor(r+0.5); % K gives the index for the row just before r
kp1 = k + 1;      % K+1 gives the index for the row just after r
r = r - k;        % R is the ratio between the K and K+1 rows

% Find indices that are out of the range 1 to n and cap them
k(k<1 | isnan(k)) = 1;
kp1 = bsxfun( @min, kp1, n );

% Use simple linear interpolation for the valid precentages
y = bsxfun(@times, 0.5-r, x(k,:)) + bsxfun(@times, 0.5+r, x(kp1,:));

% Make sure that values we hit exactly are copied rather than interpolated
exact = (r==-0.5);
if any(exact)
    y(exact,:) = x(k(exact),:);
end

function y = interpColsDiffer(x, p, n)
%INTERPCOLSDIFFER A simple 1-D linear interpolation of columns that can
%deal with columns with differing numbers of valid entries (vector n).

[nrows, ncols] = size(x);

% Make p a column vector. n is already a row vector with ncols columns.
if isrow(p)
    p = p';
end

% Form the grid of index values (numel(p) x numel(n))
r = (p/100)*n;
k = floor(r+0.5); % K gives the index for the row just before r
kp1 = k + 1;      % K+1 gives the index for the row just after r
r = r - k;        % R is the ratio between the K and K+1 rows

% Find indices that are out of the range 1 to n and cap them
k(k<1 | isnan(k)) = 1;
kp1 = bsxfun( @min, kp1, n );

% Convert K and Kp1 into linear indices
offset = nrows*(0:ncols-1);
k = bsxfun( @plus, k, offset );
kp1 = bsxfun( @plus, kp1, offset );

% Use simple linear interpolation for the valid precentages.
% Note that NaNs in r produce NaN rows.
y = (0.5-r).*x(k) + (0.5+r).*x(kp1);

% Make sure that values we hit exactly are copied rather than interpolated
exact = (r==-0.5);
if any(exact(:))
    y(exact) = x(k(exact));
end
