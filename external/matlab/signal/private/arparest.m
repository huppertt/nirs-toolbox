function [a,e,msg,msgobj] = arparest( x, p, method)
%ARPAREST   AR parameter estimation via a specified method.
%   A = ARPAREST(X,ORDER,METHOD) returns the polynomial A corresponding to 
%   the AR parametric signal model estimate of vector X using the specified
%   METHOD.  ORDER is the model order of the AR system.
%
%   Supported methods are: 'covariance' and 'modified' although all of the
%   methods of CORRMTX will work. In particular if 'autocorrelation' is
%   used, the results should be the same as those of ARYULE (but slower).
%
%   [A,E] = ARPAREST(...) returns the variance estimate E of the white noise
%   input to the AR model.

%   Ref: S. Kay, MODERN SPECTRAL ESTIMATION,
%              Prentice-Hall, 1988, Chapter 7
%        S. Marple, DIGITAL SPECTRAL ANALYSIS WITH APPLICATION,
%              Prentice-Hall, 1987, Chapter 8.
%        P. Stoica and R. Moses, INTRODUCTION TO SPECTRAL ANALYSIS,
%              Prentice-Hall, 1997, Chapter 3

%   Author(s): R. Losada and P. Pacheco
%   Copyright 1988-2014 The MathWorks, Inc.

narginchk(3,3)

% Initialize in case we return early
a = []; e = [];

% Assign msg in case there are no errors
msg ='';
msgobj = [];

% enforce backwards compatibility with row vectors
if isvector(x)
  x = x(:);
end

validateattributes(x,{'numeric'},{'nonempty','finite','2d'},'arpest','X');

% Set up necessary but not sufficient conditions for the correlation
% matrix to be nonsingular. From (Marple)
switch method,
case 'covariance',
   minlength_x = 2*p;
case 'modified',
   minlength_x = 3*p/2;
otherwise
   msgobj = message('signal:arparest:UnknMethod');
   msg = getString(msgobj);
   return
end

% Do some data sanity testing
if size(x,1) < minlength_x
   if strcmp(method, 'modified')
       multiplier = '3/2';
    else
       multiplier = '2';
    end
    if isvector(x)
       msgID = 'signal:arparest:VectorTooSmallForModel';
    else
       msgID = 'signal:arparest:MatrixTooSmallForModel';
    end
    msgobj = message(msgID,'X',multiplier);
    msg = getString(msgobj);
    return
end
if issparse(x),
   msgobj = message('signal:arparest:InputSignalCannotBeSparse');
   msg = getString(msgobj);
   return
end
if isempty(p) || p ~= round(p),
   msgobj = message('signal:arparest:ModelOrderMustBeInteger');
   msg = getString(msgobj);
   return
end

% 'like' will also copy over the complexity we only want the class
a = zeros(p+1,size(x,2),class(x)); %#ok<ZEROLIKE>
e = zeros(1,size(x,2),class(x)); %#ok<ZEROLIKE>

for chan=1:size(x,2)
   % Generate the appropriate data matrix
   XM = corrmtx(x(:,chan),p,method);
   Xc = XM(:,2:end);
   X1 = XM(:,1);

   % Coefficients estimated via the covariance method
   a(:,chan) = [1; -Xc\X1];

   % Estimate the input white noise variance
   Cz = X1'*Xc;
   e(:,chan) = X1'*X1 + Cz*a(2:end,chan);

   % Ignore the possible imaginary part due to numerical errors and force
   % the variance estimate of the white noise to be positive
   e(:,chan) = abs(real(e(:,chan)));  
end

a = a.'; % By convention all polynomials are row vectors

% [EOF] arparest.m
