function [data, truth] = simFunData( data, beta, stimLength, stimSpace, basis )
    if nargin < 5
        basis = Dictionary({'default'}, {nirs.design.basis.Canonical()});
    end
    
    if nargin < 4, stimSpace    = 10;   end
    if nargin < 3, stimLength   = 2;    end
    if nargin < 2, beta         = 5;    end
    
    if nargin < 1
        data = nirs.testing.simARNoise();
    end
    
    % sd pairs
    [SD, ~, iSD] = unique( data.probe.link(:,1:2), 'rows' );

	nSD = max(iSD);
    
    % num channels to add stim to
    nAct = floor( nSD / 2 );
    if mod(nSD, 2)
        nAct = nAct + (rand() < 0.5);
    end
    
    % choose channels
	iAct = randperm(nSD, nSD/2);
    
    % design mat
    [X, stim] = nirs.testing.randDesignMat( data.time, stimLength, stimSpace, basis );
    data.stimulus = stim;
    
    % add stimulus
    truth = zeros( size( data, 1), 1 );
    for i = 1:length(iAct)
        lst = iSD == iAct(i);
        data.data(:,lst) = data.data(:,lst) + repmat( X*beta, [1 sum(lst)] );
        
        truth = truth + lst;
    end
    
end

