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
    
    % loop through pairs
    X = []; T = table( [],[],[],[], 'VariableNames', {'source', 'detector', 'type', 'condition'} );
    for i = 1:max(iSD)
        lst = find(iSD == i);
        
        lambda = link.type(lst);

        if spectralFlag
            [xhbo, n] = nirs.design. ...
                createDesignMatrix( stimulus, t, basis, 'hbo' );

            [xhbr, n] = nirs.design. ...
                createDesignMatrix( stimulus, t, basis, 'hbr' );

            x = [];
            for j = 1:length(lambda)
                e = nirs.media.getspectra( lambda(j) );
                x = [x; [e(1)*xhbo e(2)*xhbr]];
            end
            
            src_idx = SD(i,1*ones(size(x,2),1))';
            det_idx = SD(i,2*ones(size(x,2),1))';
            mtype   = repmat({'hbo', 'hbr'}, [length(n) 1]);
            conds   = [n; n];
            
            tbl = table( src_idx, det_idx, mtype(:), conds, ...
                'VariableNames',  {'source', 'detector', 'type', 'condition'} );
        else
            error('')% test this first
            x = [];
            tbl = table( [],[],[],[], 'VariableNames', {'source', 'detector', 'type', 'condition'} );
            for j = 1:length(lambda)
                [xlam, n] = nirs.design. ...
                    createDesignMatrix( stimulus, t, basis, lambda );
                
                x = blkdiag(x, xlam);

                tbl = [tbl; 
                    table( SD(i,1*ones(size(xlam,2),1))', SD(i,2*ones(size(xlam,2),1))', ...
                    repmat(lambda(j), [size(xlam,2) 1]), ...
                    n, ...
                    'VariableNames',  {'source', 'detector', 'type', 'condition'} )];
            end
        end
        
        X   = blkdiag(X, x);
        T   = [T; tbl];
    end

    % sort
    [T, idx] = sortrows(T, {'condition', 'source', 'detector', 'type'});
    X = X(:, idx);

end

