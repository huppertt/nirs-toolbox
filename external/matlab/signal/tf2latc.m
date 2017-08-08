function [K,V]=tf2latc(varargin)
%TF2LATC Transfer function to lattice filter conversion.
%   K = TF2LATC(NUM) finds the lattice parameters K for an FIR (MA)
%   lattice filter, normalized by NUM(1).  Note that an error may
%   be generated if any zeros of the FIR filter lie on the unit circle.
%
%   K = TF2LATC(NUM,'max') where NUM corresponds to a maximum phase
%   FIR filter, reverses and conjugates NUM before converting to
%   FIR lattice.  This will result in abs(K) <= 1 and can be used
%   with LATCFILT's second output to implement the maximum phase filter.
%
%   K = TF2LATC(NUM,'min') where NUM corresponds to a minimum phase
%   FIR filter,  will result in abs(K) <= 1 and can be used
%   with LATCFILT's first output to implement the minimum phase filter.
%
%   K = TF2LATC(NUM,DEN,...) where DEN is a scalar is equivalent to
%   K = TF2LATC(NUM/DEN,...).
%
%   [K,V] = TF2LATC(NUM,DEN) finds the lattice parameters K and the ladder
%   parameters V for an IIR (ARMA) lattice-ladder filter, normalized by
%   DEN(1).  Note that an error may be generated if any poles of the
%   transfer function lie on the unit circle.
%
%   K = TF2LATC(1,DEN) finds the lattice parameters K for an IIR
%   all-pole (AR) lattice filter.  [K,V] = TF2LATC(B0,DEN) where B0 is
%   a scalar, returns a vector of ladder coefficients V.  Note that for
%   this case only the first element of V is nonzero.
%
%   EXAMPLE:
%      % Convert an all-pole IIR filter to lattice coefficients:
%      DEN = [1 13/24 5/8 1/3];
%      K = tf2latc(1,DEN);  % K will be [1/4 1/2 1/3]'
%
%   See also LATC2TF, LATCFILT, POLY2RC and RC2POLY.

% Reference:[1] S. K. Mitra, Digital Signal Processing, A Computer
%           Based Approach, McGraw-Hill, N.Y., 1998, Chapter 6.
%           [2] M. H. Hayes, Statistical Digital Signal Processing
%           and Modeling, John Wiley & Sons, N.Y., 1996, Chapter 6.
%
%   Copyright 1988-2013 The MathWorks, Inc.

narginchk(1,3);

% Parse input, num and den will be returned as column vectors
[num,den,firflag,phaseOpt] = inputparse(varargin{:});
switch firflag,
    case 1,
        [K,V] = fir2latc(num,phaseOpt,nargout);
    case 0,
        if length(num) == 1,
            [K,V] = allpole2latc(num,den);
            
        else
            % IIR filter
            [K,V] = iir2latc(num,den);
        end
end

%-------------------------------------------------------------------
function [num,den,firflag,phaseOpt] = inputparse(varargin)

% Initialize in case of early return
num = varargin{1};
den = [];   

firflag = 1;
phaseOpt = 'none';

% Default numerator
if isempty(num), num=0; end

if max(abs(num)) == 0,
    % num is zero, return
    error(message('signal:tf2latc:NeedNonZeroNum'));
end

switch nargin,
    case 2,
        if ischar(varargin{2}),
            den = 1;
            phaseOpt = varargin{2};
        else
            den = varargin{2};
        end
    case 3,
        den = varargin{2};
        phaseOpt = varargin{3};
end

% Check for valid phase options
validPhaseOpts = {'none','min','max'};
indx = find(strncmpi(phaseOpt, validPhaseOpts, length(phaseOpt)));
if isempty(indx),
    error(message('signal:tf2latc:UnknString'));
end
phaseOpt = validPhaseOpts{indx};

% Default denominator
if isempty(den), den=1; end

% Cast to enforce Precision Rules
if any([signal.internal.sigcheckfloattype(num,'single','tf2latc','NUM') ...
    signal.internal.sigcheckfloattype(den,'single','tf2latc','DEN')])
  num = single(num);
  den = single(den);
end

% Remove trailing zeros
num = removetrailzeros(num);
den = removetrailzeros(den);

% Determine if FIR
firflag = signalpolyutils('isfir',num,den);

if ~firflag && ~strcmpi(phaseOpt,'none'),
    error(message('signal:tf2latc:UnsupportedFlagFIR', 'min', 'max'));
end

% Normalize by leading den coefficient
num = num./den(1); den = den./den(1);

% Make them column vectors
num = num(:);
den = den(:);

%-------------------------------------------------------------------
function [K,V] = fir2latc(num,phaseOpt,nout)
%FIR2LATC  Convert FIR filter to lattice MA filter.

% Assign initial values in case of early return
V = [];

if nout == 2,
  V = num;
  % Cast to enforce Precision Rules
  K = zeros(length(V)-1,1,class(V)); %#ok<*ZEROLIKE>
else
    if strcmpi(phaseOpt,'max'),
        % Reverse and conjugate the polynomial to make it min. phase
        num = conj(num(end:-1:1));
    end
    
    % Now simply call poly2rc
    K = poly2rc(num);
    
    if ~strcmpi(phaseOpt,'none') && (max(abs(K)) > 1),
        warning(message('signal:tf2latc:NotMinOrMaxPhase', phaseOpt));
    end
end


%--------------------------------------------------------------------
function [K,V] = allpole2latc(num,den)

% All-pole filter, simply call poly2rc
K = poly2rc(den);
V = [num;zeros(size(K))];

%---------------------------------------------------------------------
function [K,V] = iir2latc(num,den)

% Make sure num and den are the same length:
[num,den] = eqtflength(num,den);
M = length(den);

% We still use poly2rc to compute the K's
K = poly2rc(den);
% Compute the V's recursively
% We compute the following recursion: (see Hayes, pp.306)
%                     M-1
% V(m) = num(M-m+1) - sum V(j)* conj(den_j(j-m))
%                     j=m+1
% where den_j is the denominator of jth order, the lower
% order denominators are found using the levdown function.

% We wiil use a matrix with the denominators of lower orders
% in each column, rlevinson returns this matrix
[~,tempmatrix] = rlevinson(den,1);
% Cast to enforce Precision Rules
V = zeros(M,1,class(K)); % Initialize V with zeros.
for m = M:-1:1,
    subterm = tempmatrix*V;
    V(m) = num(m) - subterm(m);
end

% EOF

