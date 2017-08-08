function DGDF = dgdfgen(Hd,hTar,doMapCoeffsToPorts)
%DGDFGEN generate the dg_dfilt structure from a specified filter structure
%Hd.

%   Author(s): Honglei Chen
%   Copyright 1988-2004 The MathWorks, Inc.


DGDF=df2tsosdggen(Hd.filterquantizer,Hd,hTar.coeffnames,doMapCoeffsToPorts,hTar.privStates);

