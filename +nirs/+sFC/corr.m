function [R,p,dfe]=corr(data,robust_flag)

if(nargin<2)
    robust_flag=true;
end

if(~isempty(strfind(class(data),'.core.Data')))
    data=data.data;
end


if(iscomplex(data))
    mask=~(imag(data)>0);
else
    mask=ones(size(data));
end


if(robust_flag)
    [R,p]=nirs.math.robust_corrcoef(data,false,mask);
else
    [R,p]=nirs.math.corrcoef(data,mask);
end

dfe = length(data);
