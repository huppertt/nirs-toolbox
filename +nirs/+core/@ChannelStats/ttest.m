function S = ttest(obj, c, b, names)
    %% ttest - Tests the null hypothesis c*beta = b
    % 
    % Args:
    %     c - contrast vector/matrix
    %     b - (optional) mean of the null distribution
    %     
    % Example:
    %     % tests the sum and difference of betas
    %     stats.ttest([1 1; 1 -1])

    nchan = size(obj.probe.link,1);

    % sort variables
    [~, icond] = sort(obj.conditions);
    obj = obj.sorted();
    
    c = c(:, icond);
    
    % full contrast matrix
    C = kron(eye(nchan), c);

    if nargin < 3
        b = zeros(size(C,1),1);
    else
        b = repmat(b(icond), [nchan 1]);
    end

    % transform beta
    beta = bsxfun(@minus, C*obj.beta, b);

    % new covariance
    covb = C*obj.covb*C';

    % output
    S = obj;

    S.beta  = beta;
    S.covb  = covb;

    % new condition names
    if nargin < 4
        cond = obj.transformNames(c);
    else
        cond = names;
    end
    
    cond = repmat( cond(:), [nchan 1] );

    link = repmat( obj.probe.link, [size(c,1) 1] );
    link = sortrows(link, {'source', 'detector', 'type'});
    S.variables = [link table(cond)];
end