function data=seperate_kernel_moments(data,keepMoment)


for i=1:length(data)
    lst=find(ismember(data(i).probe.link.moment,keepMoment));
    data(i).data=data(i).data(:,lst);
    data(i).probe.link=data(i).probe.link(lst,:);
end
