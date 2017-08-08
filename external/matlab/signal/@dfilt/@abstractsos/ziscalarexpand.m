function S = ziscalarexpand(Hd,S)
%ZISCALAREXPAND Expand empty or scalar initial conditions to a vector.

% This should be a private method

%   Author: R. Losada
%   Copyright 1988-2005 The MathWorks, Inc.

error(nargchk(2,2,nargin,'struct'));

nsecs = nsections(Hd);

% If the user passes in a cell array, it might be the old R13 style states.
if iscell(S)
    
    % Warn that this functionality might be removed.
    warning(message('signal:dfilt:abstractsos:ziscalarexpand:deprecatedFeature'));
    zicell = S;

    % Build the states from the cell array.
    S      = [];
    for indx = 1:nsecs
        S = [S zicell{indx}{2}];
    end
end
if ~isnumeric(S)
    error(message('signal:dfilt:abstractsos:ziscalarexpand:MustBeNumeric'));
end
if issparse(S),
    error(message('signal:dfilt:abstractsos:ziscalarexpand:Sparse'));
end

if nsecs ~=0
    if isempty(S),
        S = nullstate1(Hd.filterquantizer);
    end
    statespersec=2;
    if length(S)==1,
        % Zi expanded to a matrix with 2 of rows and as many columns as
        % sections
        S = S(ones(statespersec,nsecs));
    end
    % At this point we must have a matrix with the right number of rows
    if size(S,1) ~= statespersec,
        error(message('signal:dfilt:abstractsos:ziscalarexpand:InvalidDimensions', statespersec));
    end
end
