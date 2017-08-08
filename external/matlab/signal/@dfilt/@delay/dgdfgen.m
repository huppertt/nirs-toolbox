function DGDF = dgdfgen(Hd,hTar,doMapCoeffsToPorts)
%DGDFGEN Generates the dg_dfilt structure

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

DGDF = delaydggen(Hd.filterquantizer,Hd,hTar.privStates);

% [EOF]
