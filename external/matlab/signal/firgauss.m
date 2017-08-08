function [h,n] = firgauss(varargin)
%FIRGAUSS FIR Gaussian digital filter design.
%   FIRGAUSS has been replaced by GAUSSDESIGN.  FIRGAUSS still works but
%   may be removed in the future. Use GAUSSDESIGN instead. Type help
%   GAUSSDESIGN for details.
%
%   See also RCOSDESIGN.

%   Author: P. Costa
%   Copyright 1988-2013 The MathWorks, Inc.

%   References:
%     [1] Wells, W. M. III, Efficient Synthesis of Gaussian Filters by Cascaded 
%         Uniform Filters, IEEE Transactions on Pattern Analysis and Machine 
%         Intelligence, Vol. PAMI-8, No. 2, March 1980
%     [2] Rau, R., and McClellan J. H., Efficient Approximation of Gaussian Filters, 
%         IEEE Transactions on Signal Processing, Vol. 45, No. 2, February 1997.

error(nargchk(2,3,nargin,'struct'));

% Parse input arguments
[k,n] = parseinputs(varargin{:});

% Build the uniform-coefficient (boxcar) filter 
b = ones(1,n);
g = dfilt.dffir(b);

% Generate cascaded filter
Hd=cell(1,k);
[Hd{:}]=deal(g);
Hcas = cascade(Hd{:});

% Return the coefficients
h = impz(Hcas);


%-------------------------------------------------------------------
%                       Utility Functions
%-------------------------------------------------------------------
function [k,n] = parseinputs(k,varargin)
% Parse the input arguments
%
% Outputs:
%   k   - Number of cascades (always specified)
%   n   - Length of the uniform-coefficient filter
%   msg - Error message

% Default values
n = []; variance = [];
isminord   = 0;

if ischar(varargin{1}),
    isminord = 1; % N was not specified
    variance = varargin{2};
elseif isnumeric(varargin{1}),
    n = varargin{1};
end

if isminord && k >= 4,    
    % Variance was specified and number of cascades is 4 or more, 
    % compute the length, n using equation #4 in [2].
    n = getOptimalN(k,variance);
elseif isminord,
    % Variance was specified, compute the minimum n by rearranging the 
    % variance equation on page 237 of [1]. This will result in an n 
    % which meets the variance criterion specified, but is not optimal.
    n = round(sqrt((12*variance)/k + 1));
end

if k < 1
  error(message('signal:firgauss:MustBePositive', 'K'));
end
if n < 1
  error(message('signal:firgauss:MustBePositive', 'N'));
end


%-------------------------------------------------------------------
function n = getOptimalN(k,variance)
%   Determine the Optimal length for the uniform-coefficient filters
%   boxcar filters using equation 4 from [2].
%
%   Inputs:
%     k - Number of cascades  
%     variance - Overall variance of the gaussian filter
%
%   Output:
%     n - Optimal length of the uniform-coefficient filter

sigma = sqrt(variance);
if (sigma >= 2) && (sigma <= 400) && (k >= 4) && (k <= 8),
    alfa = 0.005;
else
    alfa = 0;
end

summation = 0;
for p = 0:floor(((k-1)/2)),
    summation = summation + ((-1^p/factorial(k-1))*nchoosek(k,p)*((k/2)-p)^(k-1));
end
n = abs(ceil((sqrt(2*pi)*summation + alfa)*sigma)); 

% [EOF]
