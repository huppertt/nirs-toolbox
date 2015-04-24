function stims = convertStimDesignStruct( StimDesign )

    stims = Dictionary();
    
    for i = 1:length(StimDesign)
       thisStim = convertOneStim( StimDesign(i) );
       
       % if key exists merge
       if stims.iskey( thisStim.name )
           thisStim = nirs.design.mergeStims( ...
               {stims(thisStim.name), thisStim}, ...
               thisStim.name);
       end
       
        stims(thisStim.name) = thisStim;
        
    end

end

function s = convertOneStim( stim )
    s = nirs.design.StimulusEvents();
    try
        s.name = stim.name;
    catch
        s.name = stim.cond;
    end
    
    s.onset = stim.onset(:);
    s.dur   = stim.dur(:);
    s.amp   = stim.amp(:);
end