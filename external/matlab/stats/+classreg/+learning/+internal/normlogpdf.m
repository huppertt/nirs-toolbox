function y = normlogpdf(x,mu,sigma)
    %NORMLOGPDF Normal log probability density function (pdf).
    %   Y = NORMLOGPDF(X,MU,SIGMA) returns the log pdf of the normal distribution with
    %   mean MU and standard deviation SIGMA, evaluated at the values in X.
    %   The size of Y is the common size of the input arguments.  A scalar
    %   input functions as a constant matrix of the same size as the other
    %   inputs.
    %
    %   Default values for MU and SIGMA are 0 and 1 respectively.
    %
    %   See also NORMPDF.

    %   References:
    %      [1] Evans, M., Hastings, N., and Peacock, B. (1993) Statistical
    %          Distributions, 2nd ed., Wiley, 170pp.

    %   Copyright 2014 The MathWorks, Inc. 
    
    if nargin<1
        error(message('stats:normpdf:TooFewInputs'));
    end
    if nargin < 2
        mu = 0;
    end
    if nargin < 3
        sigma = 1;
    end

    % Return NaN for out of range parameters.
    sigma(sigma <= 0) = NaN;

    try
        y = -0.5 * (((x - mu)./sigma).^2 + log(2*pi)) - log(sigma);
    catch
        error(message('stats:normpdf:InputSizeMismatch'));
    end
end
