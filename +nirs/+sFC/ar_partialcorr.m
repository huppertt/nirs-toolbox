function [R,p,dfe]=ar_partialcorr(data,modelorder,robust_flag)

if(nargin<3)
    robust_flag=false;
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


if(iscomplex(data))
    mask=~(imag(data)>0);
else
    mask=ones(size(data));
end

[yfilt,f] = nirs.math.innovations(real(data),p);

% Mask out boundary values
for ch = 1:size(yfilt,2)
    yfilt(1:length(f{ch})) = nan;
end

[i,j]=find(isnan(yfilt));
yfilt(:,unique(j))=[];

if(robust_flag)
    
    [R,p]=nirs.math.robust_partialcorr(yfilt);
else
    [R,p]=partialcorr(yfilt);
end

dfe = mean(sum(mask)) - 2;

