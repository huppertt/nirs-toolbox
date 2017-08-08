function S = ziscalarexpand(Hd,S)
%ZISCALAREXPAND Expand empty or scalar initial conditions to a vector.

% This should be a private method

%   Author: V. Pellissier
%   Copyright 1988-2005 The MathWorks, Inc.

error(nargchk(2,2,nargin,'struct'));

if ~isnumeric(S)
  error(message('signal:dfilt:singleton:ziscalarexpand:MustBeNumeric'));
end
if issparse(S),
    error(message('signal:dfilt:singleton:ziscalarexpand:Sparse'));
end

numstates = nstates(Hd);

if numstates,
    if isempty(S),
        S = nullstate1(Hd.filterquantizer);
    end
    if length(S)==1,
        % Zi expanded to a vector of length equal to the number of states
        S = S(ones(numstates,1));
    end
    
    % Transpose if row vector only.  If the filter has a single state, but
    % we have a row vector, the user probably wants to set up multichannel
    % filtering, don't transpose.
    if numstates ~= 1 && any(find(size(S)==1)),
        S = S(:);
    end
    
    % At this point we must have a vector or matrix with the right number of
    % rows
    if size(S,1) ~= numstates,
        error(message('signal:dfilt:singleton:ziscalarexpand:InvalidStateDimensions', numstates));
    end
elseif ~isempty(S),
    
    % This handles the case where one of the dimensions is zero.
    error(message('signal:dfilt:singleton:ziscalarexpand:MustBeEmpty'));
end
