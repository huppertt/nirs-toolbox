function y = unshiftdata(x,perm,nshifts)
%UNSHIFTDATA  The inverse of SHIFTDATA.
%   Y = UNSHIFTDATA(X,PERM,NSHIFTS) restores the orientation of the data that
%   was shifted with SHIFTDATA.  PERM is the permutation vector, and NSHIFTS
%   is the number of shifts that were returned from SHIFTDATA.
%
%   UNSHIFTDATA is meant to be used in tandem with SHIFTDATA.  They are handy
%   for creating functions that work along a certain dimension, like FILTER,
%   GOERTZEL, SGOLAYFILT and SOSFILT.
%
%   Examples:
%     x = magic(3)
%     [x,perm,nshifts] = shiftdata(x,2) % Work along 2nd dimension
%     y = unshiftdata(x,perm,nshifts)   % Reshapes back to original
%
%     x = 1:5                            % Originally a row
%     [x,perm,nshifts] = shiftdata(x,[]) % Work along 1st nonsingleton dimension
%     y = unshiftdata(x,perm,nshifts)    % Reshapes back to original
%
%   See also SHIFTDATA, IPERMUTE, SHIFTDIM.
  
%   Author: Thomas A. Bryan
%   Copyright 1988-2002 The MathWorks, Inc.

% Convert back to the original shape
if isempty(perm)
  y = shiftdim(x, -nshifts);
else
  y = ipermute(x,perm);
end
