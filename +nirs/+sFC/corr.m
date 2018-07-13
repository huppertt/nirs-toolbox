function [R,p,dfe]=corr(data,robust_flag)

if(nargin<2)
    robust_flag=true;
end

if(~isempty(strfind(class(data),'.core.Data')))
    data=data.data;
end


if(~isreal(data))
    mask=(imag(data)>0);
    data=real(data);
else
    mask=ones(size(data));
end


if(robust_flag)
    [R,p]=nirs.math.robust_corrcoef(data,false,mask);
else
    if(all(mask(:)==1))
        [R,p]=corrcoef(data);
    else
        [R,p]=nirs.math.corrcoef(data,false,mask);
    end
end

dfe = mean(sum(mask)) - 2;
