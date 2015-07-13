function S = ttest(obj, c, b)
    % test hypothesis that C*beta = b

    nchan = size(obj.probe.link,1);

    C = kron(eye(nchan), c);

    if nargin < 3
        b = zeros(size(C,1),1);
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
    cond = obj.transformNames(c);
    cond = repmat( cond(:)', [nchan 1] );
    cond = cond(:);

    link = repmat( obj.probe.link, [size(c,1) 1] );
    S.variables = [link table(cond)];
end