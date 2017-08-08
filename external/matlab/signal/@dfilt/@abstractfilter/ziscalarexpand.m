function S = ziscalarexpand(Hm,S)
%ZISCALAREXPAND Expand empty or scalar initial conditions to a vector.

% This should be a private method

%   Author: V. Pellissier
%   Copyright 1988-2006 The MathWorks, Inc.

error(nargchk(2,2,nargin,'struct'));

if ~isnumeric(S)
  error(message('signal:dfilt:abstractfilter:ziscalarexpand:MustBeNumeric'));
end
if issparse(S),
    error(message('signal:dfilt:abstractfilter:ziscalarexpand:Sparse'));
end

numstates = nstates(Hm);
% dss start
% if (numstates==0 && ~isempty(S)) 
%     disp('here')
%     error('signal:dfilt:abstractfilter:ziscalarexpand:NotSupported','Filter coefficients not defined. States not set.');
% end
% dss end
if numstates,
    if isempty(S),
        S=0;
    end
    if length(S)==1,
        % Zi expanded to a vector of length equal to the number of states
        S = S*ones(numstates,1);
    end
    
    % Transpose if row vector only
    if find(size(S)==1),
        S = S(:);
    end
    
    % At this point we must have a vector or matrix with the right number of
    % rows
    if size(S,1) ~= numstates,
        error(message('signal:dfilt:abstractfilter:ziscalarexpand:InvalidDimensions', numstates));
    end
elseif ~isempty(S),
    
    % This handles the case where one of the dimensions is zero.
    error(message('signal:dfilt:abstractfilter:ziscalarexpand:DFILTErr'));
end
