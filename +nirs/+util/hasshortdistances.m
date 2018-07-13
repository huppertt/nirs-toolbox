function flag = hasshortdistances(data)

for i=1:numel(data)
    flag(i)=ismember('ShortSeperation',data(i).probe.link.Properties.VariableNames);
end