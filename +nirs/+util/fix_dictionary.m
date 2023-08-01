function dictfix=fix_dictionary(dict)

dictfix=Dictionary;
st=dict.toStruct;
f=fields(st);
for i=1:length(f)
    dictfix(f{i})=st.(f{i});
end
end