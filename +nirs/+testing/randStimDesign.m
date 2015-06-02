function stim = randStimDesign( t, stimLength, stimSpace )
    
    % min max times
    tmin = min( t ) + 1*stimLength;
    tmax = max( t ) - 2*stimLength;

    % number of stims onsets
    nrnd = round( 2*(tmax-tmin)/stimSpace );
    
    % random times between tasks
    dt = exprnd(stimSpace, [nrnd 1]);

    % onsets
    onset = tmin + cumsum([0; dt]);
    onset = onset( onset < tmax );

    % durations
    dur = stimLength * ones(size(onset));

    % amplitude
    amp = ones(size(dur));

    % output
    stim = nirs.design.StimulusEvents();
    stim.amp = amp;
    stim.dur = dur;
    stim.onset = onset;
    
end

