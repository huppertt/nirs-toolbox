function Hs = sysobj(this,varargin)
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

[b, ~, ~, mssgObj] = isfdtbxinstalled;

% If an extra logical input equal to true is passed to the sysobj method,
% then sysobj returns a flag set to true if the System object conversion is
% supported for the dfilt class at hand, and to false if the conversion is
% not supported. This input is undocumented.
if nargin > 1
  validateattributes(varargin{1},{'logical'},{'scalar'},'','return flag input');
  if varargin{1}
    Hs = b && tosysobj(this,false);
    return;
  end
end

if ~b
  error(mssgObj)
end

Hs = tosysobj(this,true);

% set the meta data until the end since the System object deletes the meta
% data every time a property is set
setsysobjmetadata(this,Hs);

% [EOF]
