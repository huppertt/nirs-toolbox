function blocksetup(Hd)
%BLOCKSETUP Setup the filter object to the correct blocklength in
%preparation for calling the block method.
%
% This method is used by filterbuilder where we cannot error out and ask
% the user to change the block length of the DFILT object since the user
% does not have a handle to a DFILT object in that case. 

%   Copyright 2012 The MathWorks, Inc.

Nfft = Hd.BlockLength+Hd.ncoeffs-1;
B_Nfft = 2^nextpow2(Nfft);
if B_Nfft~=Nfft,
  Hd.BlockLength = B_Nfft-Hd.ncoeffs+1;
end
