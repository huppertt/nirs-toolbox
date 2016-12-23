function [R,p,dfe]=partialcorr(data,robust_flag)

if(nargin<2)
    robust_flag=true;
end

if(~isempty(strfind(class(data),'.core.Data')))
    data=data.data;
end

if(iscomplex(data))
    warning('code does not support masked data yet');
    data=real(data);
end

if(robust_flag)
    warning('robust partial corr not yet supported')
    robust_flag=false;
end


if(robust_flag)
    [R,p]=nirs.math.robust_corrcoef(data);
else
    [R,p]=partialcorr(data);
end

dfe = length(data);
