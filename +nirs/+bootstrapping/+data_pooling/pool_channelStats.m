function result = pool_channelStats(data,type)

if(nargin<2)
    type='beta';
end

for i=1:length(data)
    if(strcmp(type,'beta'))
        result(:,i)=data(i).beta;
    elseif(strcmp(type,'tvalue'))
        result(:,i)=data(i).tvalue;
    else
        result(:,i)=data(i).(type);
    end
end


