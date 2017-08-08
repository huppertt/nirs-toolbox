%SYSOBJ Generate a filter System object
%   Hs = sysobj(Hd) generates a filter System object, Hs, based on the
%   dfilt object, Hd.
%
%   The supported structures are:
%    
%   FIR Structures (sysobj generates a dsp.FIRFilter System object):
%     Direct form FIR (dfilt.dffir)
%     Direct form FIR transposed (dfilt.dffirt)
%     Direct form symmetric FIR (dfilt.dfsymfir)
%     Direct form anti symmetric FIR (dfilt.dfasymfir)
%     Lattice MA for min phase (dfilt.latticemamin)
%
%   SOS Structures (sysobj generates a dsp.BiquadFilter System object):
%     Direct form I (dfilt.df1sos)
%     Direct form I transposed (dfilt.df1tsos)
%     Direct form II (dfilt.df2sos)
%     Direct form II transposed (dfilt.df2tsos)


%   Copyright 2011 The MathWorks, Inc.

% [EOF]
