function [X, stim] = randDesignMat( t, stimLength, stimSpace, basis )
    if nargin < 5
        basis = Dictionary({'default'}, {nirs.design.basis.Canonical()});
    end

    stim = nirs.testing.randStimDesign( t, stimLength, stimSpace );

    stim = Dictionary( {'rand'}, {stim} );
    
    % design mat
    X = nirs.design.createDesignMatrix( stim, t, basis );
end