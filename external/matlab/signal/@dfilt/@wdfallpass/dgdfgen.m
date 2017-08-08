function DGDF = dgdfgen(Hd,hTar,doMapCoeffsToPorts)
%DGDFGEN generate the dg_dfilt structure from a specified filter structure
%Hd.

%   Author(s): Honglei Chen
%   Copyright 1988-2005 The MathWorks, Inc.


DGDF=wdfallpassdggen(Hd.filterquantizer,Hd,hTar.coeffnames,doMapCoeffsToPorts,hTar.privStates);

