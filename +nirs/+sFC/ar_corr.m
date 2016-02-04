function [R,p,dfe]=ar_corr(data,modelorder,robust_flag)

if(nargin<3)
    robust_flag=true;
end

if(nargin<2 || isempty(modelorder))
    modelorder=4;
end

if(isa(data,'nirs.core.Data'))
    Fs=data.Fs;
    data=data.data;
else
    Fs=1;
end


if(isstr(modelorder))
    p = Fs*str2num(modelorder(1:strfind(modelorder,'x')-1));
else
    p=modelorder;
end

[yfilt,f] = nirs.math.innovations(data,p);

if(robust_flag)
    [R,p]=nirs.math.robust_corrcoef(yfilt);
else
    [R,p]=corrcoef(yfilt);
end

dfe = length(yfilt);
