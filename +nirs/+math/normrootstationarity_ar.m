function dataf = normrootstationarity_ar(data)


md=mean(data,1);
data=data-ones(size(data,1),1)*md;
vd=var(data,[],1);
data=data./(ones(size(data,1),1)*sqrt(vd));

[innov,f]=nirs.math.innovations(data,10);

innov_n=nirs.math.normrootstationarity(innov,'std');
innov_n=innov_n./(ones(size(innov_n,1),1)*sqrt(var(innov_n,[],1)./var(innov,[],1)));

dataf=data;
for i=1:size(data,2)
    dataf(:,i)=filter(1,f{i},innov_n(:,i));
end

dataf=dataf-ones(size(dataf,1),1)*mean(dataf,1);
dataf=dataf./(ones(size(dataf,1),1)*sqrt(var(dataf,[],1)));
dataf=dataf.*(ones(size(dataf,1),1)*sqrt(vd));
dataf=dataf+ones(size(dataf,1),1)*md;


