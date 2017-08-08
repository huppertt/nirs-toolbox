function [Ht,anum,aden] = iirbpc2bpc(Hd, varargin)
%IIRBPC2BPC IIR complex bandpass to complex bandpass transformation

%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.

[Ht,anum,aden] = ciirxform(Hd, @zpkbpc2bpc, varargin{:});

