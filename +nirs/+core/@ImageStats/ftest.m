function S = ftest(obj, m)
    %% ftest - This uses Hotelling's T squared test to calculate F-stats
    % 
    % Args:
    %     m - each row specifies a mask of the conditions that 
    %         should be jointly tested
    %         
    % Example:
    %     % performs two differents F-tests of all 3 variables and also
    %     % of the first two variables
    %     stats.ftest([1 1 1; 1 1 0])
    
    if ~islogical(m)
        m = m > 0;
        warning('Converting mask to true/false.')
    end

    nchan = size(obj.beta,1)/length(obj.conditions);

    % sort variables
    [~, icond] = sort(obj.conditions);
    obj = obj.sorted();
    
    m = m(:, icond);
    
    % full contrast matrix
    M = kron(speye(nchan,nchan), m);

    n = obj.dfe * ones(size(M,1),1);
    k = sum(M,2);

    df1 = k;
    df2 = n-k+1; % adding 1 so that df2 matches tstats dfe

    for i = 1:size(M,1)
        idx = M(i,:) > 0;
        b = obj.beta(idx);
        
        covb=obj.covb_chol(idx,:)*obj.covb_chol(idx,:)';
        
        T2(i,1) = b'*pinv(covb)*b;
    end

    F = (n-k+1) ./ k ./ n .* T2;

    % new condition names
    cond = obj.transformNames(m);
    for i = 1:length(cond)
        cond{i} = strjoin(strsplit(cond{i},'+'),' & ');
    end
    cond = repmat( cond(:), [nchan 1] );

%     link = repmat( obj.probe.link, [size(m,1) 1] );
%     link = sortrows(link, {'source', 'detector', 'type'});
    
    S = nirs.core.ImageFStats();
    S.variables = obj.variables; 
   
    S.F         = F;
    S.df1       = full(df1);
    S.df2       = full(df2);
    S.probe     = obj.probe;
    S.mesh      = obj.mesh;
end