function stim = blockedStimDesign( t, stimDur, stimSpace, ncond )
%%blockedStimDesign - Returns an alternating blocked design task.
% 
% Args:
%     t           - time vector
%     stimDur     - stim duration
%     stimSpace   - average space between stim onsets
%     ncond       - number of conditions

    if nargin < 4
        ncond = 1;
    end
    
    onset = rand(1,1)*20;  % random onset to blocks
    
    % figure out num-trials
    timepertrial = (stimDur+stimSpace)*ncond;
    numtrials = floor((max(t)-min(t)-onset)/timepertrial);
    
    % output
    stim = Dictionary();
    
    onsetsall=(min(t)+onset):(stimDur+stimSpace):max(t)-stimDur;
    
    for i = 1:ncond
        s = nirs.design.StimulusEvents();
        s.onset = onset + onsetsall(i:ncond:end);
        s.dur   = stimDur*ones(size(s.onset));
        s.amp   = ones(size(s.onset));        
        s.name   = char(uint8('A') + i - 1);
        
        stim(s.name) = s;
    end
    
end

