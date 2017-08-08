function c = sos2cell(s,g)
%SOS2CELL Convert second-order-section matrix to cell array.
%   C = SOS2CELL(S) converts L-by-6 second-order-section matrix S 
%   in the form
%
%        S =   [B1 A1
%               B2 A2
%                ...
%               BL AL]
%
%   to a cell array C in the form
%
%     C = { {B1,A1}, {B2,A2}, ... {BL,AL} }
%
%   where each numerator vector Bi and denominator vector Ai represents the
%   coefficients of a linear or quadratic polynomial. 
%
%   C = SOS2CELL(S,G) with an additional gain term prepends a constant term to C
%   in the form
%
%     C = { {G,1}, {B1,A1}, {B2,A2}, ... {BL,AL} }
%
%   Example:
%     [b,a] = butter(4,.5);
%     [s,g] = tf2sos(b,a);
%     c = sos2cell(s,g)
%
%   See also CELL2SOS, TF2SOS, SOS2TF, ZP2SOS, SOS2ZP, SOS2SS, SS2SOS.

%   Thomas A. Bryan
%   Copyright 1988-2004 The MathWorks, Inc.

if nargin<2
  g = [];
end

if ~isempty(g) && g == 1;
    % Assume there is no g if the gain is 1
    g = [];
end

if ~isnumeric(s) || ndims(s)~=2 || size(s,2)~=6,
   error(message('signal:sos2cell:InvalidDimensions'))
end
m = size(s, 1);

% If we have a single-section with [1 0 0 1 0 0], recombine the gain
if m == 1 && ~isempty(g) && ~any((s==[1 0 0 1 0 0]) ~= 1) ,
    s(1) = g*s(1);
    g = [];
end

isg = ~isempty(g);
c = cell(1,m+isg); % m for S plus 1 for g if not empty
if isg
  c{1} = {g,1};
end
for i=1:m
  % Assign cells after removing trailing zeros
  c{i+isg}={dezero(s(i,1:3)), dezero(s(i,4:6))};
end
 
%%%%%%%%%%%%%%%%%%%%%%%%
function x1 = dezero(x)
%DEZERO Remove trailing zeros.
%   DEZERO(X) removes trailing zeros from numeric vector X.
%

%   Modified from DEBLANK

if ~isempty(x) && ~isnumeric(x) && ~isvector(x)
    warning(message('signal:sos2cell:InvalidParam'))
end

if isempty(x)
    x1 = x([]);
else
  % remove trailing zeros
  [r,c] = find(x ~= 0); %#ok
  if isempty(c)
    x1 = x([]);
  else
    x1 = x(1:max(c));
  end
end
 
