function onset = first_event(data)

onset=zeros(length(data),1);
for i=1:length(data)
    on=[];
    for s=1:data(i).stimulus.count;
        ss=data(i).stimulus(data(i).stimulus.keys{s});
        if(isa(ss,'nirs.design.StimulusEvents'))
            on=[on; ss.onset(:)];
        end
    end
    onset(i)=min(on);
end
    