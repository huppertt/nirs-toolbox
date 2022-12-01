function Dnew=FixDictionaryCount(D)

Dnew=Dictionary;
for i=1:D.count; 
    Dnew(D.keys{i})=D.values{i}; 
end;