function stims = convertStimDesignStruct( StimDesign,time )

if(nargin<2)
    time=[];
end



    stims = Dictionary();
    
    for i = 1:length(StimDesign)
       thisStim = convertOneStim( StimDesign(i),time );
       
       % if key exists merge
       if stims.iskey( thisStim.name )
           thisStim = nirs.design.mergeStims( ...
               {stims(thisStim.name), thisStim}, ...
               thisStim.name);
       end
       
        stims(thisStim.name) = thisStim;
        
    end

end

function s = convertOneStim( stim,time )

if(~isfield(stim,'continuous') || ~stim.continuous)
    
    
    s = nirs.design.StimulusEvents();
    try
        s.name = stim.name;
    catch
        s.name = stim.cond;
    end
    
    s.onset = stim.onset(:);
    
    if length(stim.dur) == 1
        s.dur = stim.dur * ones(size(stim.onset));
    else
        s.dur   = stim.dur(:);
    end
    
    if(~isfield(stim,'amp'))
        stim.amp=1;
    end
    
    if length(stim.amp) == 1
        s.amp = stim.amp * ones(size(stim.onset));
    else
        s.amp = stim.amp(:);
    end
else
    s = nirs.design.StimulusVector();
       try
        s.name = stim.name;
    catch
        s.name = stim.cond;
    end
    s.vector=stim.amp(:);
    s.time=time(:);
    
end
end