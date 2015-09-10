function S = ttest(obj, c, b)
    %% ttest - Tests the null hypothesis c*beta = b
    % 
    % Args:
    %     c - contrast vector/matrix
    %     b - (optional) mean of the null distribution
    %     
    % Example:
    %     % tests the sum and difference of betas
    %     stats.ttest([1 1; 1 -1])

    nchan = length(obj.beta)/length(obj.conditions);

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
    covb=C*obj.covb_chol;
  
    % output
    S = obj;

    S.beta  = beta;
    S.covb_chol  = covb;
    S.typeII_StdE = C*obj.typeII_StdE;
    
    % new condition names
    cond = obj.transformNames(c);
    cond = repmat( cond(:), [nchan 1] );

    var = obj.variables;
    type=unique(var.cond);
    lst=find(ismember(var.cond,type{1}));
    var=var(lst,:);
    var.cond=cond;
    S.variables = var;
end