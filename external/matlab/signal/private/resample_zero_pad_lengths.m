function [nz,nz1,delay] = resample_zero_pad_lengths(p,q,Lx,L)
%RESAMPLE_ZERO_PAD_LENGTHS computes the delays for the RESAMPLE function
%
%   [NZ,NZ1,DELAY] = RESAMPLE_ZERO_PAD_LENGTHS(P,Q,LX,L) returns the number
%   of leading zeros (NZ), trailing zeros (NZ1) that the RESAMPLE
%   function uses to pad its input to UPFIRDN, and the amount of DELAY that
%   must be stripped from the output.
%
%   This is a private function used by toolbox/signal/eml/resample.m.

% Copyright 2009 The MathWorks, Inc.

    Lhalf = (L-1)/2;
    nz = floor(q-mod(Lhalf,q));
    
    delay = floor(ceil(Lhalf+nz)/q);

    % Need to zero-pad so output length is exactly ceil(Lx*p/q).
    nz1 = 0;
    while ceil( ((Lx-1)*p+L+nz1 )/q ) - delay < ceil(Lx*p/q)
        nz1 = nz1+1;
    end

    
