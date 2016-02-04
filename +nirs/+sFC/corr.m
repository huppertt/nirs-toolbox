function [R,p,dfe]=corr(data,robust_flag)

if(nargin<2)
    robust_flag=true;
end

if(isa(data,'nirs.core.Data'))
    data=data.data;
end

if(robust_flag)
    [R,p]=nirs.math.robust_corrcoef(data);
else
    [R,p]=corrcoef(data);
end

dfe = length(data);
