function [X, T] = createMultiVarDesignMat( stimulus, t, basis, link, spectralFlag )
    % assert probe is sorted
    [~,idx] = sortrows(link, {'source', 'detector', 'type'});
    
    if ~all( diff(idx) == 1 )
       error( 'Please sort channels by source/detector/type.')
    end

    if nargin < 5
        spectralFlag = true;
    end
    
    % unique sd pairs
	[SD, ~, iSD] = unique([link.source link.detector], 'rows');
    
    % loop through source-detector pairs
    X = []; 
    T = table( [],[],[],[], 'VariableNames', {'source', 'detector', 'type', 'cond'} );
    
    for i = 1:max(iSD)
        
        % list of channels for pair i
        lst = find(iSD == i);
        
        % wavelengths
        lambda = link.type(lst);
        
        if spectralFlag && ~isnumeric(lambda)
            error('Multivariate design with spectral priors requires optical density')
        end

        % with spectral priors?
        if spectralFlag
            
            % get the oxy-hb design matrix
            [xhbo, n] = nirs.design. ...
                createDesignMatrix( stimulus, t, basis, 'hbo' );

            % get the deoxy design matrix
            [xhbr, n] = nirs.design. ...
                createDesignMatrix( stimulus, t, basis, 'hbr' );

            % weight the columns by extinction coef
            x = [];
            for j = 1:length(lambda)
                e = nirs.media.getspectra( lambda(j) );
                x = [x; [e(1)*xhbo e(2)*xhbr]];
            end
            
            % need a table to keep track of all the columns of X
            src_idx = SD(i,1*ones(size(x,2),1))';
            det_idx = SD(i,2*ones(size(x,2),1))';
            mtype   = repmat({'hbo', 'hbr'}, [length(n) 1]);
            conds   = [n; n];
            
            tbl = table( src_idx, det_idx, mtype(:), conds, ...
                'VariableNames',  {'source', 'detector', 'type', 'cond'} );
            
        % no spectral priors (dOD)
        else
            
            % loop through each wavelength
            x = [];
            tbl = table( [],[],[],[], 'VariableNames', {'source', 'detector', 'type', 'cond'} );
            for j = 1:length(lambda)
                
                % get design matrix for lambda j
                [xlam, n] = nirs.design. ...
                    createDesignMatrix( stimulus, t, basis, lambda );
                
                x = blkdiag(x, xlam);
                
                % table to keep track of columns of X
                src_idx = SD(i,1*ones(size(xlam,2),1))';
                det_idx = SD(i,2*ones(size(xlam,2),1))';
                
                mtype = repmat(lambda(j), [size(xlam,2) 1]);
                
                tbl = [tbl; 
                    table( src_idx, det_idx, mtype, n, ...
                    'VariableNames',  {'source', 'detector', 'type', 'cond'} )];
            end
        end
        
        % add the contribution of source-detector pair i to
        % the design matrix
        X   = blkdiag(X, x);
        T   = [T; tbl];
    end

    % sort using the convention i chose
    [T, idx] = sortrows(T, {'source', 'detector', 'type', 'cond'});
    X = X(:, idx);

end

