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

    nchan = size(obj.probe.link,1);

    M = kron(m, eye(nchan)) > 0;

    n = obj.dfe * ones(size(M,1),1);
    k = sum(M,2);

    df1 = k;
    df2 = n-k;

    for i = 1:size(M,1)
        idx = M(i,:);
        b = obj.beta(idx);
        T2(i,1) = b'*pinv(obj.covb(idx,idx))*b;
    end

    F = (n-k) ./ k ./ (n-1) .* T2;


    % new condition names
    condition = obj.transformNames(m);
    for i = 1:length(condition)
        condition{i} = strjoin(strsplit(condition{i},'+'),' & ');
    end

    condition = repmat( condition(:)', [nchan 1] );
    condition = condition(:);

    link = repmat( obj.probe.link, [size(m,1) 1] );
    vars = [link table(condition)];

    S = nirs.core.ChannelFStats();

    S.variables = vars;
    S.F         = F;
    S.df1       = df1;
    S.df2       = df2;
    S.probe     = obj.probe;
end