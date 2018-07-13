function replace_stimdata(data)

for i=1:length(data)
    StimDesign=struct('name',[],'onset',[],'dur',[],'amp',[]);
    for j=1:data(i).stimulus.count
        s=data(i).stimulus(data(i).stimulus.keys{j});
         StimDesign(j).name=s.name;
         StimDesign(j).onset=s.onset;
         StimDesign(j).dur=s.dur;
         StimDesign(j).amp=s.amp;
    end
    tmp=load(data(i).description,'-MAT');
    if(isfield(tmp,'StimDesign'))
        tmp.StimDesignOld=tmp.StimDesign;
    end
    tmp.StimDesign=StimDesign;
    save(data(i).description,'-STRUCT','tmp','-MAT');
    disp(data(i).description)
end