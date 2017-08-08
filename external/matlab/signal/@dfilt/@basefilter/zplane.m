function [z,p,k] = zplane(Hb,varargin)
%ZPLANE Z-plane zero-pole plot.
%   ZPLANE(Hb) plots the zeros and poles of the discrete-time filter Hb with the
%   unit circle for reference.  Each zero is represented with a 'o' and each
%   pole with a 'x' on the plot.  Multiple zeros and poles are indicated by the
%   multiplicity number shown to the upper right of the zero or pole.
%
%   See also DFILT.   
  
%   Copyright 1999-2004 The MathWorks, Inc.

Hd = dispatch(Hb);
if nargout>0,
    [z,p,k] = zpk(Hd);
else
    [Hb, opts] = freqzparse(Hb, varargin{:});
    fvtool(Hb, 'polezero', opts);
end

% [EOF]
