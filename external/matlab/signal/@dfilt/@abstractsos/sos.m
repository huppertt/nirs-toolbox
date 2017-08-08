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
%   See also DFILT.

%   Author: R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

error(nargchk(1,2,nargin,'struct')); % No scaling allowed

if nargin < 2,
    Hsos = copy(Hd);
else
    % Convert to second-order-section matrix and gain.
    [z,p,k] = zpk(Hd);
    [s,g] = zp2sos(z,p,k,varargin{:});
    
    % Use same structure as this Hd.
    Hsos = feval(str2func(class(Hd)),s,g);
end
