function lphpreorder(this,Hd)
%LPHPREORDER   Rule-of-thumb lowpass/highpass reordering of SOS.

%   Author(s): R. Losada
%   Copyright 1999-2004 The MathWorks, Inc.

% Reference: D. Schlichtharle. Digital Filters. Basic and Design.
% Springer-Verlag. Berlin, 2000.

reorder(Hd,'down');

nsecs = nsections(Hd);

% Get reorder indices for lowpass/highpass rule-of-thumb
reorderindx = lphpreorderindx(this,Hd,nsecs);

reorder(Hd,reorderindx);



% [EOF]
