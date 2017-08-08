function [grid,Ngrid,edge_vec,edge_idx] = crmz_grid(edge_pairs, L, fdomain, ...
                                                grid_density)
% CRMZ_GRID Generate grid for Remez exchange.
%
%   [Grid, Ngrid, V_edges, Indx_edges] = ...
%                    CRMZ_GRID(Edge_Pairs, L, fdomain, grid_density)
%           Grid: frequencies on the fine grid
%          Ngrid: number of grid pts if whole domain used [-1,+1]
%       edge_vec: specified edges reshaped into a vector of adjacent edge-pairs
%       edge_idx: index of specified band-edges within grid array
%
%     edge_pairs: specified band-edges, [Nbands X 2] array
%              L: filter length
%        fdomain: domain for frequency approximation
%                 default is [0,0.5] called 'half'
%                 'whole' means [-0.5,0.5), 
%   grid_density: density of grid

%   Authors: L. Karam, J. McClellan
%   Copyright 2005 The MathWorks, Inc.

error(nargchk(3,4,nargin,'struct'));

if (edge_pairs(1) == -0.5) && (edge_pairs(end) == 0.5),
   % -pi and +pi are the same point - move one a little:
   new_freq = 0.5 - 1/(50*L);
   if new_freq <= edge_pairs(end-1),
     % Last two freq points are too close to move - try first two:
     new_freq = -new_freq;
     if (new_freq >= edge_pairs(2)),
       error(message('signal:crmz_grid:InvalidParam'));
     else
       edge_pairs(1) = new_freq;
     end
   else
     edge_pairs(end) = new_freq;
   end
end

Ngrid = 2.^nextpow2(L*grid_density);
if (Ngrid/2) > 20*L,
   Ngrid = Ngrid/2;
end

edge_vec = edge_pairs';  % M-by-2 to 2-by-M
edge_vec = edge_vec(:);  % single column of adjacent edge-pairs

switch fdomain
case 'whole'
  grid = (0:Ngrid-1)./Ngrid - 0.5;             % uniform grid points [-.5,.5)
  edge_idx = 1 + round((edge_vec+0.5)*Ngrid);  % closest indices in grid
case 'half'
  grid = (0:Ngrid/2)./Ngrid;                   % uniform grid points [0,.5]
  edge_idx = 1 + round(edge_vec*Ngrid);        % closest indices in grid
otherwise
  error(message('signal:crmz_grid:InternalError'));
end
edge_idx(end) = min(length(grid), edge_idx(end));  % Clip last index

% Fix repeated edges:
% This determines not only if
%         i(1)==i(2) and i(3)==i(4),
% but also if
%         i(2)==i(3), etc.
m = find(edge_idx(1:end-1) == edge_idx(2:end));
if ~isempty(m),
  % Replace REPEATED band edges with the uniform grid points
  % Could be a problem if [-1 -1] (if whole) or [0 0] (if half) specified
  edge_idx(m) = edge_idx(m) - 1;   % move 1 index lower
  edge_vec(m) = grid(edge_idx(m)); % change user's band edge accordingly
  m = m + 1;
  edge_idx(m) = edge_idx(m) + 1;
  edge_vec(m) = grid(edge_idx(m));
end

% Replace closest grid points with exact band edges:
grid(edge_idx) = edge_vec;

% end of crmz_grid

