function dataf = normrootstationarity(data,Stat)
% Stat
%      'mean'   detect changes in mean
%      'rms'    detect changes in root-mean-square level
%      'std'    detect changes in standard deviation
%      'linear' detect changes in mean and slope


if(nargin<2 || isempty(Stat))
    Stat='std';
end
dataf=data;
for i=1:size(data,2)
    [a]=findchangepts(data(:,i),'Statistic',Stat);
    
    a=[1; a; size(data,1)];
    for j=1:length(a)-1
        dataf(a(j):a(j+1)-1,i)=renorm(dataf(a(j):a(j+1)-1,i),Stat);
    end

end


function d = renorm(d,Stat)

switch(Stat)
    case('mean')
        d=d-median(d);
    case('std')
        md=median(d);
        d=(d-md)./mad(d,1)+md;
    case('rms')
        d=d-sqrt(median(d.^2));
    case('linear')
        x=[[1:length(d)]' ones(length(d),1)];
        b=nirs.math.robustfit(x,d,[],[],'off');
        d=d-x*b;
end
