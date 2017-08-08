function Hsos = sos(Hd,varargin)
%SOS  Convert to second-order-sections.
%   Hsos = SOS(Hd) converts discrete-time filter Hd to second-order section
%   form.
%
%   SOS(Hd,DIR_FLAG) specifies the ordering of the 2nd order sections. If
%   DIR_FLAG is equal to 'UP', the first row will contain the poles closest
%   to the origin, and the last row will contain the poles closest to the
%   unit circle. If DIR_FLAG is equal to 'DOWN', the sections are ordered
%   in the opposite direction. The zeros are always paired with the poles
%   closest to them. DIR_FLAG defaults to 'UP'.
%
%   SOS(Hd,DIR_FLAG,SCALE) specifies the desired scaling of the gain
%   and the numerator coefficients of all 2nd order sections. SCALE can be
%   either 'NONE', Inf or 2 which correspond to no scaling, infinity
%   norm scaling and 2-norm scaling respectively. SCALE defaults to 'NONE'.
%   The filter must be stable in order to scale in the 2-norm or inf-norm sense.
%   Using infinity-norm scaling in conjunction with 'UP' ordering will
%   minimize the probability of overflow in the realization. On the other
%   hand, using 2-norm scaling in conjunction with 'DOWN' ordering will
%   minimize the peak roundoff noise.
%
%   See also DFILT.

%   Author: R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

error(nargchk(1,3,nargin,'struct')); % Scaling allowed

if nargin < 2,
    Hsos = copy(Hd);
else
    % Convert to second-order-section matrix and gain.
    [z,p,k] = zpk(Hd);
    [s,g] = zp2sos(z,p,k,varargin{:});
    
    % Use same structure as this Hd.
    Hsos = copy(Hd);
    Hsos.sosmatrix = s;
    Hsos.ScaleValues = g;
end
