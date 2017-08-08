function bpbsreorder(this,Hd)
%BPBSREORDER  Rule-of-thumb bandpass/bandstop reordering of SOS.

%   Author(s): R. Losada
%   Copyright 1999-2004 The MathWorks, Inc.

% Reference: D. Schlichtharle. Digital Filters. Basic and Design.
% Springer-Verlag. Berlin, 2000. 

reorder(Hd,'down');

nsecs = nsections(Hd);

% Get reorder indices for lowpass/highpass rule-of-thumb
lpreorderindx = lphpreorderindx(this,Hd,ceil(nsecs/2));

% Initialize reorderindx vector
reorderindx = zeros(1,2*size(lpreorderindx,2));

% Set odd sections
reorderindx(1:2:end) = 2*lpreorderindx-1;

% Set even sections
reorderindx(2:2:end) = 2*lpreorderindx;

% Eliminate extra section if odd nsections
if rem(nsecs,2) == 1,
    mindx = find(reorderindx == max(reorderindx));
    reorderindx(mindx) = [];
end

reorder(Hd,reorderindx);


% [EOF]
