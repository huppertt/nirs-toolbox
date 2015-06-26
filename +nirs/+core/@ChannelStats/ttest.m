function S = ttest(obj, c, b)
    % test hypothesis that C*beta = b

    nchan = size(obj.probe.link,1);

    C = kron(c, eye(nchan));

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
    condition = obj.transformNames(c);
    condition = repmat( condition(:)', [nchan 1] );
    condition = condition(:);

    link = repmat( obj.probe.link, [size(c,1) 1] );
    S.variables = [link table(condition)];
end