function demo = combine_demographics(demoOrig)

demo=Dictionary;
flds = demoOrig.Properties.VariableNames;
for i=1:length(flds)
    if(length(unique(demoOrig.(flds{i})))==1)
        demo(demoOrig.Properties.VariableNames{i})=demoOrig.(flds{i})(1);
    end
end



end