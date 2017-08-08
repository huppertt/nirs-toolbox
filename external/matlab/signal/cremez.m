function [h,delta,result] = cremez(M, edges, filt_str, varargin)
%CREMEZ Complex and nonlinear phase equiripple FIR filter design.
%   CREMEZ is obsolete.  CREMEZ still works but may be removed in the future.
%   Use CFIRPM instead.
%
%   See also CFIRPM.

%   Authors: L. Karam, J. McClellan
%   Revised: October 1996, D. Orofino
%
%   Copyright 1988-2004 The MathWorks, Inc.

% NOTE: This algorithm is equivalent to Remez for real B
%       when the filter specs are exactly linear phase.

if nargin<3, error(message('signal:cremez:Nargchk')); end

[h,delta,result] = cfirpm(M, edges, filt_str, varargin{:});

% [EOF]
