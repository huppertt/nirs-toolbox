function [A,B,C,D] = ss(Hd) %#ok<STOUT,INUSD>
%SS  Discrete-time filter to state-space conversion (not supported).
%   dfilt.allpass does not support the analysis method ss. 
%   For other dfilt objects, ss works as follows:
%   [A,B,C,D] = SS(Hd) converts discrete-time filter Hd to state-space
%   representation given by 
%     x(k+1) = A*x(k) + B*u(k)
%     y(k)   = C*x(k) + D*u(k)
%   where x is the state vector, u is the input vector, and y is the output
%   vector. 

%   Copyright 2013-2014 The MathWorks, Inc.

% Return error message 'Method ss not supported by dfilt.allpass filters.'
error(message('signal:dfilt:allpass:allpass:noStateSpaceSupport'))
