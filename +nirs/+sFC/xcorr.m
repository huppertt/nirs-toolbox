function [R,p,dfe]=xcorr(data,maxlags,robust_flag)

if(nargin<2)
    maxlags=[];
end

if(nargin<3)
    robust_flag=true;
end

if(~isempty(strfind(class(data),'.core.Data')))
    Fs=data.Fs;
    data=data.data;
else
    Fs=1;
end

if(ischar(maxlags))
    maxlags = ceil(Fs*str2double(maxlags(1:strfind(maxlags,'x')-1)));
end

if(~isreal(data))
    mask=(imag(data)>0);
    data=real(data);
else
    mask=ones(size(data));
end

if robust_flag
    [R,p]=nirs.math.robust_xcorrcoef(data,maxlags,mask);
else
    [R,p]=nirs.math.xcorrcoef(data,maxlags,mask);
end

dfe = mean(sum(mask)) - 2;
