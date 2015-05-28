function [X, stim] = randDesignMat( t, stimLength, stimSpace, basis )
    if nargin < 5
        basis = Dictionary({'default'}, {nirs.design.basis.Canonical()});
    end

    % random stim design
    tmin = min( t ) + stimLength;
    tmax = max( t ) - 2*stimLength;

    nrnd = round( 2*(tmax-tmin)/stimSpace );
    dt = exprnd(stimLength, [nrnd 1]);

    onset = tmin + cumsum([0; dt]);
    onset = onset( onset < tmax );

    dur = stimLength * ones(size(onset));

    amp = ones(size(dur));

    stim = nirs.design.StimulusEvents();
    stim.amp = amp;
    stim.dur = dur;
    stim.onset = onset;

    stim = Dictionary({'rand'},{stim});

    % add to data
    X = nirs.design.createDesignMatrix( stim, t, basis );
end

% X = nirs.testing.randDesignMat( linspace(0,300,1200)', 5, 10 );
% 
% Y = [1*X (-1*X)];
% 
% e = mvnrnd([0 0], [1 0.5; 0.5 1], 1200);
% 
% e(:,1) = filter(1, [1 -0.5], e(:,1));
% e(:,2) = filter(1, [1 -0.3], e(:,2));
