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

if(~isreal(data))
    mask=(imag(data)>0);
    
else
    mask=ones(size(data));
end

[yfilt,f] = nirs.math.innovations(real(data),p);


if(robust_flag)
    [R,p]=nirs.math.robust_corrcoef(yfilt,false,mask);
else
    [R,p]=nirs.math.corrcoef(yfilt,mask);
end

dfe = length(yfilt);

