function ic = super_getinitialconditions(Hd)
%SUPER_GETINITIALCONDITIONS Get the initial conditions

%   Copyright 2009 The MathWorks, Inc.

s  = Hd.States;

nic = double(s.Numerator);
dic = double(s.Denominator);

nstages = size(Hd.sosMatrix,1);         % number of SOS stages
nst = 2;                                % number of num's or den's states
% per one SOS stage
nchannels = numel(nic)/(nst*nstages);   % number of channels

ICNum = zeros(nst*nstages,nchannels);
ICDen = zeros(nst*nstages,nchannels);
for c=1:nchannels
    
    % get NIC and DIC for each channel
    startidx = ((c-1)*nstages)+1;
    stopidx = startidx+nstages-1;
    nic_channel = nic(:,startidx:stopidx);
    dic_channel = dic(:,startidx:stopidx);
    
    ICNum(:,c) = reshape(nic_channel, nst*nstages, size(nic, 1)/2);
    ICDen(:,c) = reshape(dic_channel, nst*nstages, size(dic, 1)/2);
end

ic.Num = ICNum;
ic.Den = ICDen;

% [EOF]