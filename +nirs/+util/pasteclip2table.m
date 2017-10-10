function tbl = pasteclip2table

data=clipboard('paste');

header=textscan(data(1:min(find(double(data)==13))),'%s','Delimiter','\t');
C={}; cnt=1;
s=[];
for i=1:length(header{1})
    s=[s '%s'];
end
a=textscan(data,s,'Delimiter','\t');

if(strcmp(a{1}(2,end),'-'))
    n=3;
else
    n=2;
end
s=struct;
for i=1:length(header{1})
    
    if(~isempty(str2num(a{i}{n})))
        s=setfield(s,header{1}{i},cell2mat(arrayfun(@(x)str2num(x{1}),vertcat(a{i}(n:end)),'UniformOutput',false)));
    else
        s=setfield(s,header{1}{i},vertcat(a{i}(n:end)));
    end
end


tbl=struct2table(s);

return