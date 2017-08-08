function varargout = signalpolyutils(varargin)
%SIGNALPOLYUTILS   utility functions for vectors of polynomial coefficients.
%   S = SIGNALPOLYUTILS provides access to a number of local functions that
%   have a variety of polynomial manipulation and testing utilities.

%   Author(s): R. Losada
%   Copyright 1988-2012 The MathWorks, Inc.

[varargout{1:max(1,nargout)}] = feval(varargin{:});

%----------------------------------------------------------------------------
function polynom = imagprune(polynom,tol) %#ok
%IMAGPRUNE  Remove imaginary part when smaller than tol.

if nargin<2, tol=[]; end
if isempty(tol), tol = eps^(2/3); end
 
if max(abs(imag(polynom))) < tol,
   polynom = real(polynom);
end

%-----------------------------------------------------------------------------
function symstr = symmetrytest(b,removezerosFlag,tol)
%SYMMETRYTEST  Test if vector corresponds to a symmetric or antisymmetric polynomial.

if nargin<2, removezerosFlag=[]; end
if isempty(removezerosFlag), removezerosFlag = 0; end

if nargin<3, tol=[]; end
if isempty(tol), tol = eps^(2/3); end

% Make sure b is a row
b = b(:).';

if removezerosFlag,
    % Remove leading and trailing zeros of b 
    b = removezeros(b);
end

% Try complex first
switch isreal(b)
case 0,
    if max(abs(b - conj(b(end:-1:1)))) <= tol,
        symstr = 'hermitian';
    elseif max(abs(b + conj(b(end:-1:1)))) <= tol,
        symstr = 'antihermitian';
    else
        symstr = 'none';
    end
    
case 1,
    if max(abs(b - b(end:-1:1))) <= tol,
        symstr = 'symmetric';
    elseif max(abs(b + b(end:-1:1))) <= tol,
        symstr = 'antisymmetric';
    else
        symstr = 'none';
    end    
end

%------------------------------------------------------------------------------------
function filtertype = determinetype(h,issymflag,removezerosFlag) %#ok
%DETERMINETYPE  Determine the type of the filter based on 
%               the length and the symmetry of the filter.

if removezerosFlag,
    % Remove leading and trailing zeros of b 
    h = removezeros(h);
end
N = length(h) - 1;

if issymflag,
    % Type 1 or type 2
    if rem(N,2),
        % Odd order
        filtertype = 2;
    else
        % Even order
        filtertype = 1;
    end
else
    % Type 3 or type 4
    if rem(N,2),
        % Odd order
        filtertype = 4;
    else
        % Even order
        filtertype = 3;
    end
end

%-----------------------------------------------------------------------------
function flag = isminphase(b,tol)
%ISMINPHASE  Test to see if polynomial has all its roots on or inside the unit circle.

if nargin<2, tol=[]; end
if isempty(tol), tol = eps^(2/3); end

flag = 1;

% First test if polynomial is strictly minimum-phase, i.e. all its roots are
% strictly inside the unit circle.
stableflag = isstable(b);
if stableflag,
    return
end

% If not strictly minimum-phase, it can still be minimum-phase, try this.

% Remove trailing zeros of b before calling roots, otherwise, the order of
% the input polynomial will be incorrect. 
if ~isempty(b)
  b1 = b(1:find(b~=0, 1, 'last'));
else
  b1 = b;
end

z = roots(b1);
if ~isempty(z) && (max(abs(z)) > 1 + tol),
    flag = 0;        
end

%------------------------------------------------------------------------------
function flag = isstable(a)
%ISSTABLE  Test to see if polynomial has all its roots inside the unit circle.

% Remove trailing zeros
a = a(1:find(a~=0, 1, 'last' ));

% Remove leading zeros as they have no effect on stability but affect the
% normalization
indx = find(a, 1);
if ~isempty(indx),
    a = a(indx:end);
else
    % All zeros
    error(message('signal:signalpolyutils:SignalErr'));
end

a = a./a(1);    % Normalize by a(1)

if length(a) == 1,
    flag = 1;
    
elseif length(a) == 2,
    
    % One pole given by second coefficient of a, first is always 1.
    flag = isfirstorderstable(a(2));

else
    % Use poly2rc for denominators of order 2 or more
    flag = ispolystable(a);
        
end

%----------------------------------------------------------------------------
function isstableflag = isfirstorderstable(p)
% One pole given by second coefficient of a, first is always 1.
if abs(p) < 1,        
    isstableflag = 1;
else
    isstableflag = 0;
end

%----------------------------------------------------------------------------
function isstableflag = ispolystable(a)

% look at the last coefficient, if greater or equal to 1, unstable
if abs(a(end)) >= 1,
    isstableflag = 0;
    return
end

% Use poly2rc to determine stability
try
    k = poly2rc(a); % This can throw a divide by zero warning
    if any(isnan(k)) || max(abs(k)) >= 1,
        isstableflag = 0;
    else
        isstableflag = 1;
    end
catch %#ok<CTCH>
    % If poly2rc fails, one of the k's must be equal to one, unstable
    isstableflag = 0;
end


%----------------------------------------------------------------------------
function isfirflag = isfir(b,a)
%ISFIR(B,A) True if FIR.
if nargin<2, a=[]; end
if isempty(a), a=1; end
if ~isvector(b) || ~isvector(a)
  error(message('signal:signalpolyutils:InvalidDimensions'));
end

if find(a ~= 0, 1, 'last') > 1,
  isfirflag = 0;
else
  isfirflag = 1;
end

%----------------------------------------------------------------------------
function islinphaseflag = islinphase(b,a,tol) %#ok
%ISLINPHASE(B,A) True if linear phase
if nargin<3, tol=[]; end
if isempty(tol), tol=eps^(2/3); end
if isfir(b,a),
  islinphaseflag = determineiflinphase(b,tol);
else
  if isstable(a),
    % Causal stable IIR filters cannot have linear phase
    islinphaseflag = 0;
  else
    islinphaseflag = (determineiflinphase(b,tol) & determineiflinphase(a,tol));
  end
end        

%----------------------------------------------------------------------------
function islinphaseflag = determineiflinphase(b,tol)
if nargin<2, tol=[]; end
if isempty(tol), tol=eps^(2/3); end

% Set defaults
islinphaseflag = 0;

% If b is a scalar, filter is always FIR and linear-phase
if length(b) == 1,
    islinphaseflag = 1;
    return
end

symstr = symmetrytest(b,1,tol);

if ~strcmpi(symstr,'none'),
    islinphaseflag = 1;
end
 
%----------------------------------------------------------------------------
function ismaxphaseflag = ismaxphase(b,a,tol) %#ok
%ISMAXPHASE(B,A) True if maximum phase
% Initialize flag to true.
if nargin<3, tol=[]; end
if isempty(tol), tol = eps^(2/3); end
ismaxphaseflag = 1;

% Remove trailing zeros of b before calling roots, otherwise, the order of
% the input polynomial will be incorrect. 
if ~isempty(b)
  b1 = b(1:find(b~=0, 1, 'last'));
else
  b1 = b;
end

% If there is a zero at the origin, filter is not max phase
if any(roots(b1)==0) || length(a) > length(b) || ...
     ~isminphase(conj(b(end:-1:1)),tol) || ...
     ~isstable(a);
  ismaxphaseflag = 0;
end

%----------------------------------------------------------------------------
function isallpassflag = isallpass(b,a,tol) %#ok
%ISALLPASS(B,A,TOL) returns true if allpass.
% If the numerator and denominator are conjugate reverses of each other,
% then the filter is allpass. 

%isallpass(b,a,tol)  True if allpass.
if nargin<3, tol=[]; end
if isempty(tol), tol = eps^(2/3); end

% Get out if empty so no indexing errors occur.
if isempty(a) && isempty(b)
  isallpassflag = 1;
  return
end

% Remove trailing zeros.
b = removetrailzeros(b);
a = removetrailzeros(a);

% Remove leading zeros.
b = b(find(b, 1):end);
a = a(find(a, 1):end);

% Get out if one of them is now empty so no indexing errors occur.
if isempty(a) || isempty(b)
  isallpassflag = 0;
  return
end

% Get out of one is zero so no divide-by-zero warnings occur.
if a(1)==0 || b(end)==0
  isallpassflag = 0;
  return
end

% Normalize
a = a./a(1);
b = b./b(end);

% If the numerator and denominator are conjugate reverses of each other,
% then the filter is allpass. 
if length(b)==length(a) && max(abs(conj(b(end:-1:1)) - a)) < tol,
  isallpassflag = 1;
else
  isallpassflag = 0;
end

%----------------------------------------------------------------------------
function t = isvector(v)
%ISVECTOR  True for a vector.
%   ISVECTOR(V) returns 1 if V is a vector and 0 otherwise.
t = ismatrix(v) & min(size(v))<=1;

%----------------------------------------------------------------------------
function b = removezeros(b)
%REMOVEZEROS  Remove leading and trailing zeros of b

if max(abs(b)) == 0,
    b = 0;
else
    % Remove leading and trailing zeros of b 
    b = b(find(b,1):find(b,1, 'last'));
end
