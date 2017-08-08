function Q = marcumq(varargin)
%MARCUMQ Generalized Marcum Q function.
%   MARCUMQ(A,B,M) is the generalized Marcum Q function, defined as
%
%      Q_m(a,b) = 1/a^(m-1) * integral from b to inf of
%                 [x^m * exp(-(x^2+a^2)/2) * I_(m-1)(ax)] dx,
%
%   where I_(m-1)() is the modified Bessel function of the first kind, of
%   order m-1.
%
%   a and b must be real and nonnegative.  m must be a positive integer. 
%   If any of the inputs is a scalar, it is expanded to the size of the 
%   other inputs.
%
%   MARCUMQ(A,B) is the special case for M=1, originally tabulated by
%   Marcum and often written without the subscript, Q(a,b). 
%
%   The calculation in the function is based on the method proposed by
%   Shnidman, using an absolute error criterion.
%
%   % Example:
%   %   Compute Marcum Q function value of 10 and 20.
%   
%   y = marcumq(10,20)
%
%   See also BESSELI.

%   Copyright 1996-2012 The MathWorks, Inc.

%   References:
%     [1] D. A. Shnidman, "The Calculation of the Probability of Detection
%         and the Generalized Marcum Q-Function", IEEE Transactions on
%         Information Theory, vol. IT-35, pp. 389-400, March 1989.
%     [2] J. I. Marcum, "A Statistical Theory of Target Detection
%         by Pulsed Radar:  Mathematical Appendix", RAND Corporation, Santa
%         Monica, CA, Research Memorandum RM-753, 1 July 1948. Reprinted in
%         IRE Transactions on Information Theory, vol. IT-6, pp. 59-267,
%         April 1960.

% Check number of input arguments
narginchk(2,3);
% Assign input arguments to variables a, b, and m
if nargin==2, a = varargin{1}; b = varargin{2}; m=1; end
if nargin==3, a = varargin{1}; b = varargin{2}; m=varargin{3}; end

% Check the input data type. Single precision is not supported.
try
    chkinputdatatype(a,b,m);
catch ME
    throwAsCaller(ME);
end

if (any(~isreal(a)) || any(~isreal(b)) || any(~isreal(m)))
    error(message('signal:marcumq:nonrealArgs'));
end

if (any(a(:)<0) || any(b(:)<0))
    error(message('signal:marcumq:negativeArgs'));
end

if (length(a)==1)
    a = repmat(a,size(b));
elseif (length(b)==1)
    b = repmat(b,size(a));
elseif (~isequal(size(a),size(b)))
    error(message('signal:marcumq:inconsistentDims'))
end

if (length(m)==1)
    m = repmat(m,size(a));
else
    if (length(a)==1)
        a = repmat(a,size(m));
        b = repmat(b,size(m));
    elseif (~isequal(size(a),size(m)))
        error(message('signal:marcumq:inconsistentDims'))
    end
end

if ~isempty(m)
    if (any(m(:)<1) || any(m(:)./floor(m(:))~=1))
        error(message('signal:marcumq:improperM'))
    end
end

% Switch to Shnidman's notation[1].  Note that because the values of a, b,
% and m will affect the choice of power series form and the limits used in
% the summations evaluated in the P function, we process the vector or
% matrix one element at a time

Q = NaN(size(a));
for idxx = 1:numel(a),
    % Special cases
    if a(idxx)~=Inf && b(idxx)==0, Q(idxx) = 1; end
    if m==1
        if a(idxx)~=Inf && b(idxx)==Inf, Q(idxx) = 0; end
        if a(idxx)==Inf && b(idxx)~=Inf, Q(idxx) = 1; end
    end
    % General case if the value is still NaN
    % Convert a and b to X and Y using Eq. 7 of [1]
    if isnan(Q(idxx)),
        N = m(idxx);
        X = (a(idxx)^2/2)/N;
        Y = b(idxx)^2/2;
        Q(idxx) = P(N,X,Y);
    end
end

%-------------------------------------------------------------------------------
% Subfunctions
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
function Q = P(N,X,Y)
% Computes the probability of detection as given by Shnidman[1].

NX = N*X;  % Notational convenience

% Evaluate the Chernoff bound
if Y==0
    logCB = 0;  % Avoid division by zero
else
    Lambda = 1-N/(2*Y)-sqrt((N/(2*Y))^2+NX/Y);  % Eq. 22
    logCB = -Lambda*Y+NX*Lambda/(1-Lambda)-N*log(1-Lambda);  % Eq. 23
end

if exp(logCB) < realmin
    if Y < NX + N
        Q = 1;
        return;
    elseif Y > NX + N
        Q = 0;
        return;
    end
else
    if Y < NX + N   % Compute 1-P(N,X,Y) using form 4 
        Form = 4;   % Form 4 corresponds to Eq. 11
    else            % Compute P(N,X,Y) using form 2 
        Form = 2;   % Form 2 corresponds to Eq. 9
    end
    
    switch Form
        case 2
            % Form 2 corresponds to Eq. 9
            % Adjust start of inner summation to avoid underflow
            kmin = findStartOfSummation(NX, 0);

            % Adjust start of outer summation to avoid underflow
            mmin = findStartOfSummation(Y, kmin+N);

            % Initialize m and k to adjusted values
            m = mmin;
            k = m-N;

            % Initialize the values of innerSum and outerSum
            %
            % If value of k based on N and m is greater than the minimum value
            % computed above, pre-compute start of inner summation.  If not,
            % then we can ignore the first k-1 terms in the inner summation.
            if kmin<k
                innerSum = 0;
                for idx = kmin:k
                    innerArg = expA(NX,idx);
                    if innerArg > realmin
                        innerSum = innerSum + innerArg;
                    end
                end
            else
                innerArg = expA(NX,k);
                innerSum = innerArg;
            end
            % With m = mmin, ignore the first m-1 terms in the outer summation
            outerArg = expA(Y,m);
            outerSum = outerArg*(1-innerSum);

            % Sum until the outerSum is no longer increasing
            while innerArg>eps(innerSum) && outerArg*(1-innerSum)>realmin
                m = m + 1;
                k = k + 1;
                innerArg = innerArg*NX/k;
                innerSum = innerSum + innerArg;
                outerArg = outerArg*Y/m;
                outerSum = outerSum + outerArg*(1-innerSum);
            end

            % Form 2 also has a summation over m separate from the double sum term
            num_million = floor(N/1e6);
            for m = 1:num_million
                outerSum = outerSum + sum(expA(Y,(m-1)*1e6:m*1e6-1));
            end
            outerSum = outerSum + sum(expA(Y,num_million*1e6:N-1));

            Q = outerSum;

        case 4
            % Form 4 corresponds to Eq. 11
            % Adjust start of inner summation to avoid underflow
            mmin = findStartOfSummation(Y, 0);

            % Adjust start of outer summation to avoid underflow
            kmin = findStartOfSummation(NX, max(mmin-N+1,0));

            % Initialize k and m to adjusted values
            k = kmin;
            m = N - 1 + k;

            % Initialize the values of innerSum and outerSum
            %
            % If value of m based on N and k is greater than the minimum value
            % computed above, pre-compute start of inner summation.  If not,
            % then we can ignore the first m-1 terms in the inner summation.
            if mmin<m
                innerSum = 0;
                for idx = mmin:m
                    innerArg = expA(Y,idx);
                    if innerArg > realmin
                        innerSum = innerSum + innerArg;
                    end
                end
            else
                innerArg = expA(Y,m);
                innerSum = innerArg;
            end
            % With k = kmin, ignore the first k-1 terms in the outer summation
            outerArg = expA(NX,k);
            outerSum = outerArg*(1-innerSum);

            % Sum until the outerSum is no longer increasing
            while innerArg>eps(innerSum) && outerArg*(1-innerSum)>realmin
                m = m + 1;
                k = k + 1;
                innerArg = innerArg*Y/m;
                innerSum = innerSum + innerArg;
                outerArg = outerArg*NX/k;
                outerSum = outerSum + outerArg*(1-innerSum);
            end

            Q = 1-outerSum;
    end
end
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
function startIdx = findStartOfSummation(constTerm, minVal)

epsM = 1e-40;

G = -4 * log(4 * epsM * (1-epsM)) / 5;  % Eq. 41

%the parameter constTerm should be large enough for Shnidman's assumptions
%to hold. If the parameter is not large enough, the ensuing calculations
%are meaningless (we get a square root with negative input in Eq. 40).
%Shnidman only offers a solution to the underflow region problem for large
%constTerm because: 
%1. This is effectively the case of interest 
%2. If constTerm is not large, the solution can't be obtained in a expedient
% manner, and is thus inefficient. 
% Here, we assume that if (2*constTerm - G)>0, the assumptions hold and the
% underflow region is investigated. Else, this step is skipped.
if (expA(constTerm,minVal) > realmin) || ( (2*constTerm - G) < 0)
    % If there is no underflow in the minVal term, or if simplifying
    % assumptions do not hold, start summation at minVal
    startIdx = minVal;
else
    % Underflow is detected in the minVal term, so compute new starting index
    startIdx = floor(constTerm + 1/2 - sqrt(G * (2*constTerm - G))); % Eq. 40
end
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
function x = expA(y,n)
% Evaluates terms of the form exp(A) = (e^-y)*(y^n)/(n!).  For large values of y
% we use a modified expression for exp(A) from a follow-up paper by Shnidman:
%
% [3] D. A. Shnidman, "Note on 'The Calculation of the Probability of
% Detection and the Generalized Marcum Q-Function'", IEEE Transactions
% on Information Theory, vol. 37, no. 4, p. 1233, July 1991.

if y>0
    if y>1e4
        % Note:  There is a typo in the equation for "A" in [3], affecting the 
        % first term: (z+0.5) should read (z-0.5) 
        % and the term in the denominator of the first term in the brackets:
        % 1+1/(2z) should read 1-1/(2z).
        z = n + 1;
        x = exp((z-0.5).*((1-y./z)./(1-1./(2*z))+log(y./z))-0.5*log(2*pi*y)-J(z));
    else
        x = exp(-y + n*log(y) - gammaln(n+1));
    end
else
    % In this case y equals zero, so the answer is zero for all cases except
    % when n equals zero, in which case the answer is one
    if n == 0
        x = 1;
    else
        x = 0;
    end
end
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
function x = J(z)
% A continued fraction approximation to the Binet function given in [1]

x = 1./(12*z+2./(5*z+53./(42*z+1170./(53*z+53./z))));   % Eq. B5
%--------------------------------------------------------------------------
%-----