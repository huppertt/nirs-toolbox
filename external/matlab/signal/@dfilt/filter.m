%FILTER Execute ("run") discrete-time filter.
%   Y = FILTER(Hd,X) filters the data X using the discrete-time filter
%   object Hd to create the filtered data Y. The final conditions are
%   stored in Hd.States. 
%
%   If Hd.PersistentMemory is false (default), initial conditions
%   are set to zero before filtering. 
%
%   To use non-zero initial conditions, set Hd.PersistentMemory to true
%   and set Hd.States  with a vector of NSTATES(Hd) elements. If a scalar
%   is specified, it will be expanded to a vector of the correct length. 
%
%   If X is a matrix, each column is filtered as an independent channel. In
%   this case, initial conditions can be optionally specified for each
%   channel individually by setting Hd.States to a matrix of NSTATES(Hd)
%   rows and SIZE(X,2) columns.  If X is a multidimensional array, filter
%   operates on the first nonsingleton dimension.
%
%   FILTER(Hd,X,DIM) operates along the dimension DIM. If X is a vector or 
%   matrix and DIM is 1, every column of X is treated as a channel. If DIM
%   is 2, every row represents a channel.
%  
%   EXAMPLE:
%     x = randn(100,1);       % Original signal
%     b = fir1(50,.4);        % 50th-order linear-phase FIR filter
%     Hd = dfilt.dffir(b);    % Direct-form FIR implementation
% 
%     % No initial conditions
%     y1 = filter(Hd,x);      % PersistentMemory is false (default)
%     zf = Hd.States;         % Final conditions
%   
%     % Non-zero initial conditions
%     Hd.PersistentMemory = true;
%     Hd.States = 1;          % Uses scalar expansion
%     y2 = filter(Hd,x);
%     stem([y1 y2])           % Different sequences at the beginning
% 
%     % Streaming data
%     reset(Hd);              % Clear filter history
%     y3 = filter(Hd,x);      % Filter the entire signal in one block
%     reset(Hd);              % Clear filter history
%     yloop = [];
%     xblock = reshape(x,[20 5]);
%     % Filtering the signal section by section is equivalent to filtering the
%     % entire signal at once.
%     for i=1:5,
%       yloop = [yloop; filter(Hd,xblock(:,i))];
%     end
% 
%   See also DFILT/NSTATES

%   Author: V. Pellissier
%   Copyright 1988-2006 The MathWorks, Inc.

% Help for the DFILT filter method.

% [EOF]
