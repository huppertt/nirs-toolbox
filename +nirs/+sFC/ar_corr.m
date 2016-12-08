function [R,p,dfe]=ar_corr(data,modelorder,robust_flag)

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


if(isstr(modelorder))
    p = Fs*str2num(modelorder(1:strfind(modelorder,'x')-1));
else
    p=modelorder;
end


lst=all(data==0,2);

[yfilt,f] = nirs.math.innovations(data,p);

yfilt(lst,:)=[];

if(robust_flag)
    [R,p]=nirs.math.robust_corrcoef(yfilt);
else
    [R,p]=corrcoef(yfilt);
end

dfe = length(yfilt);

