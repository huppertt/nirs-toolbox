function val = dictionaryGet(dict,names,default)

lst=find(ismember(lower(dict.keys),names));
if(isempty(lst))
    val=default;
else
    val=dict(dict.keys{lst(1)});
end