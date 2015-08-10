function stim = randStimDesign( t, stimDur, stimSpace, ncond )
%% randStimDesign - Returns a random stim design.
% 
% Args:
%     t           - time vector
%     stimDur     - stim duration
%     stimSpace   - average space between stim onsets
%     ncond       - number of conditions

    if nargin < 4
        ncond = 1;
    end
    
    % min max times
    tmin = min( t ) + 1*stimDur;
    tmax = max( t ) - 2*stimDur;

    % number of stims onsets
    nrnd = round( 2*(tmax-tmin)/stimSpace );
    
    % random times between tasks
    dt = stimSpace/2 + exprnd(stimSpace/2, [nrnd 1]);

    % onsets
    onset = tmin + cumsum([0; dt]);
    onset = onset( onset < tmax );

    % durations
    dur = stimDur * ones(size(onset));

    % amplitude
    amp = ones(size(dur));

    % output
    stim = Dictionary();
    r = randi(ncond);
    
    for i = 1:ncond
        lst = mod(i+r,ncond)+1:ncond:size(dur,1);
        
        s = nirs.design.StimulusEvents();
        s.name   = char(uint8('A') + i - 1);
        s.amp    = amp(lst);
        s.dur    = dur(lst);
        s.onset  = onset(lst);
        
        stim(s.name) = s;
    end
    
end

