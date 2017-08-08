function s = norm(Hd,pnorm,tol)
%NORM   Filter norm.
%   NORM(H) returns the L2-norm of a digital filter (DFILT) H.  
%
%   NORM(H,PNORM) returns the p-norm of a filter. PNORM can be either
%   frequency-domain norms: 'L1', 'L2', 'Linf' or discrete-time-domain
%   norms: 'l1', 'l2', 'linf'. Note that the L2-norm of a filter is equal
%   to its l2-norm (Parseval's theorem), but this is not true for other
%   norms.
%
%   When computing the l1-, l2-, linf-, L1-, and L2-norms of an IIR filter,
%   NORM(...,TOL) will specify the tolerance for greater or less accuracy.
%   By default, TOL = 1e-8.
%
%      EXAMPLE(S):
%      % Compute the L2-norm with a tolerance of 1e-10 for an IIR filter
%      Hs = fdesign.lowpass; % Create a filter design specifications object
%      Hd = design(Hs,'butter'); % Design a Butterworth SOS filter
%      L2 = norm(Hd,'L2',1e-10); % Compute the L2-norm
%
%   See also DFILT/SCALE, DFILT/SCALECHECK, NORM.

%   Author(s): R. Losada
%   Copyright 2003-2005 The MathWorks, Inc.

