function out = mergeStims( stims, name )
    
    onset   = [];
    dur     = [];
    amp     = [];
    
    % loop through stims and concat onset, dur, amp
    for i = 1:length(stims)
        onset   = [onset; stims{i}.onset(:)];
        dur     = [dur; stims{i}.dur(:)];
        amp     = [amp; stims{i}.amp(:)];
    end
    
    % create output
    out = nirs.functional.StimulusEvents();
    
    out.name    = name;
    out.onset   = onset;
    out.dur     = dur;
    out.amp     = amp;
end