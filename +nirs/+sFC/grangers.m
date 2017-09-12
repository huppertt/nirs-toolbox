function [G,F,dfe1,dfe2,p]=grangers(data,modelorder,robust_flag,includeZeroLag)
% This function runs the Granger's causality models
if(nargin<4)
    includeZeroLag=true;
end

if(nargin<3)
    robust_flag=true;
end

if(nargin<2 || isempty(modelorder))
    modelorder=20;
end

if(~isempty(strfind(class(data),'.core.Data')))
    Fs=data.Fs;
    data=data.data;
else
    Fs=1;
end


if(ischar(modelorder))
    Pmax = ceil(Fs*str2double(modelorder(1:strfind(modelorder,'x')-1)));
else
    Pmax = modelorder;
end

if robust_flag
    [G, F, dfe1, dfe2, p] = nirs.math.robust_mvgc(data,Pmax,includeZeroLag);
else
    [G, F, dfe1, dfe2, p] = nirs.math.mvgc(data,Pmax,includeZeroLag);
end

end