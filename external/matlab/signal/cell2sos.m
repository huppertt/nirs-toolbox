function [s,g] = cell2sos(c)
%CELL2SOS Convert cell array to second-order-section matrix
%   S = CELL2SOS(C) converts cell array C in the form
%
%     C = { {B1,A1}, {B2,A2}, ... {BL,AL} },
%
%   where each numerator vector Bi and denominator vector Ai represents the
%   coefficients of a linear or quadratic polynomial, to an L-by-6
%   second-order-section matrix S of the form
%
%        S = [B1 A1
%             B2 A2
%              ...
%             BL AL]
%
%   Linear sections are zero-padded on the right.
%
%   [S,G] = CELL2SOS(C) when the first element of C is a pair of
%   scalars, will form the scalar gain G with that pair and use the
%   remaining elements of C to build the S matrix.
%
%   Example:
%     % Gain is embedded in the leading first-order section:
%     c = {{[0.0181 0.0181],[1.0000 -0.5095]},{[1 2 1],[1 -1.2505  0.5457]}};
%     s = cell2sos(c)
%
%     % Gain is embedded in the leading zeroth-order (scalar) section:
%     c = {{0.0181,1},{[1 1],[1.0000 -0.5095]},{[1 2 1],[1 -1.2505  0.5457]}};
%     [s,g] = cell2sos(c)
%
%   See also SOS2CELL, TF2SOS, SOS2TF, ZP2SOS, SOS2ZP, SOS2SS, SS2SOS.

%   Thomas A. Bryan
%   Copyright 1988-2004 The MathWorks, Inc.

% Check the input data type. Single precision is not supported.
try
    for n = 1:length(c)
        chkinputdatatype(c{n}{1},c{n}{2});
    end
catch ME
    throwAsCaller(ME);
end

if nargout == 2,
    c1 = c{1};
    if length(c1{1}) == 1 && length(c1{2}) == 1,
        g = c1{1}./c1{2};
        
        % Use the remaining part of c to build s
        c = {c{2:end}};
    else
        g = 1;
    end
    
    
end

s = formsos(c);

%--------------------------------------------------------------
function s = formsos(c)
m = length(c);
s = zeros(m,6);
for i=1:m
  b = c{i}{1};
  a = c{i}{2};
  b=b(:).';
  a = a(:).';
  b=[b,zeros(1,3-length(b))];
  a=[a,zeros(1,3-length(a))];
  s(i,:) = [b a];
end
