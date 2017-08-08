function [n_out, w, trivalwin] = check_order(n_in)
%CHECK_ORDER Checks the order passed to the window functions.
% [N,W,TRIVALWIN] = CHECK_ORDER(N_ESTIMATE) will round N_ESTIMATE to the
% nearest integer if it is not already an integer. In special cases (N is
% [], 0, or 1), TRIVALWIN will be set to flag that W has been modified.

%   Copyright 1988-2002 The MathWorks, Inc.

w = [];
trivalwin = 0;

if ~(isnumeric(n_in) & isfinite(n_in)),
    error(message('signal:check_order:InvalidOrderFinite', 'N'));
end

% Special case of negative orders:
if n_in < 0,
   error(message('signal:check_order:InvalidOrderNegative'));
end

% Check if order is already an integer or empty
% If not, round to nearest integer.
if isempty(n_in) | n_in == floor(n_in),
   n_out = n_in;
else
   n_out = round(n_in);
   warning(message('signal:check_order:InvalidOrderRounding'));
end

% Special cases:
if isempty(n_out) | n_out == 0,
   w = zeros(0,1);               % Empty matrix: 0-by-1
   trivalwin = 1; 
elseif n_out == 1,
   w = 1;
   trivalwin = 1;   
end
