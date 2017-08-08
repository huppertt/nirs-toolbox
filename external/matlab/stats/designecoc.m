function M = designecoc(K,dname,varargin)
%DESIGNECOC Coding matrix for reducing a multiclass problem to a set of binary problems.
%   M=DESIGNECOC(K,NAME) returns a K-by-L matrix M for K classes and L
%   binary classifiers. The number of binary classifiers is set by K and
%   the value of the coding matrix design NAME. Pass K as a positive
%   integer and NAME as a string. DESIGNECOC forms M using -1, +1, and 0.
%   If M(I,J) is -1, class I is treated as negative for classifier J. If
%   M(I,J) is +1, class I is treated as positive for classifier J. If
%   M(I,J) is 0, class I is not used for training classifier J. You can
%   pass this matrix to FITCECOC to fit a multiclass model by
%   error-correcting output code (ECOC).
%
%   Allowed values of NAME and the respective numbers of classifiers for K
%   classes are summarized below:
%       -------------------------------------------------------------------
%       |       Name                |           L           |   Note      |
%       |---------------------------|-----------------------|-------------|
%       | 'onevsone' or 'allpairs'  |      K*(K-1)/2        |             |
%       | 'onevsall'                |           K           |             |
%       | 'binarycomplete'          |      2^(K-1)-1        |             |
%       | 'ternarycomplete'         | (3^K - 2^(K+1) + 1)/2 |             |
%       | 'ordinal'                 |          K-1          |             |
%       | 'sparserandom'            |      15*log2(K)       | Approximate |
%       | 'denserandom'             |      10*log2(K)       | Approximate |
%       -------------------------------------------------------------------
%
%   M=DESIGNECOC(K,NAME,'PARAM1',val1,'PARAM2',val2,...) specifies optional
%   parameter name/value pairs:
%       'NumTrials'     - Positive integer indicating the number of trials
%                         to be used for generation of the coding matrix
%                         for 'sparserandom' and 'denserandom' designs.
%                         DESIGNECOC generates that many matrices and
%                         selects the one with maximal pairwise row
%                         distance. If you pass NAME other than
%                         'sparserandom' or 'denserandom', DESIGNECOC
%                         ignores this parameter. Default: 10000
%
%   See also fitcecoc.

%   Copyright 2014 The MathWorks, Inc.

if ~isnumeric(K) || ~isscalar(K) || round(K)~=K || K<1
    error(message('stats:designecoc:BadK'));
end

if K==1
    M = 1;
    return;
end

if K==2
    M = [-1 1]';
    return;
end

allowedVals = {'onevsone' 'allpairs' 'onevsall' 'binarycomplete' 'ternarycomplete' ...
    'ordinal' 'sparserandom' 'denserandom'};

tf = strncmpi(dname,allowedVals,length(dname));
Nfound = sum(tf);
if Nfound~=1
    error(message('stats:designecoc:BadDesignName',sprintf(' ''%s''',allowedVals{:})));
end
dname = allowedVals{tf};
if strcmp(dname,'allpairs')
    dname = 'onevsone';
end

tf = ismember({'onevsone' 'onevsall' 'binarycomplete' 'ternarycomplete' ...
    'ordinal' 'sparserandom' 'denserandom'},dname);

args = {'numtrials'};
defs = {        1e4};
N = internal.stats.parseArgs(args,defs,varargin{:});

if ~isscalar(N) || ~isnumeric(N) || N<1 || N~=round(N)
    error(message('stats:designecoc:BadNumTrials'));
end

if     tf(1) % one vs one
    L = K*(K-1)/2;
    M = zeros(K,L);
    
    l = 1;
    t = 1;
    r = K-1;
    b = K;
    
    for k=1:K-1
        M(t,l:r)     = +1;
        M(t+1:b,l:r) = -eye(b-t);
        
        l = r+1;
        t = t+1;
        r = r+K-k-1;
        b = t+K-k-1;
    end
    
elseif tf(2) % one vs all
    M = -ones(K);
    M(1:K+1:end) = +1;
    
elseif tf(3) % binary complete: 2^(k-1)-1 columns
    M = ff2n(K-1);
    M(1,:) = [];
    M = [zeros(size(M,1),1) M];
    M(M==0) = -1;
    M = -M';
    
elseif tf(4) % ternary complete: (3^k - 2^(k+1) + 1)/2 columns
    M = fullfact(repmat(3,1,K))';
    M(M==1) = -1;
    M(M==2) =  0;
    M(M==3) =  1;
    M = clean(M);
    
elseif tf(5) % ordinal
    M = -ones(K,K-1);
    for k=1:K-1
        M(k+1:end,k) = 1;
    end
    
elseif tf(6) % sparse random
    if K<5
        warning(message('stats:designecoc:NotEnoughClassesForSparseRandom'));
    end

    % Ternary Hamming distance
    distfun = @(y,M) nansum( 1 - bsxfun(@times,M,y), 2 )/2;
    
    M = [];
    maxmindist = -Inf;
    
    % Set the number of binary learners (columns)
    [~,L] = log2(K);
    L = 15*L;
    
    % Generate many matrices and choose the one maximizing the minimal
    % Hamming distance between rows. Elements -1 and +1 are generated with
    % probability 1/4 each. Element 0 is generated with probability 1/2.
    for n=1:N
        Mtry = rand(K,L);
        iminus = Mtry<1/4;
        iplus  = Mtry>3/4;
        izero = ~(iminus | iplus);
        Mtry(iminus) = -1;
        Mtry(iplus)  = +1;
        Mtry(izero)  = NaN;

        % Compute pairwise distances for all rows and store in vector D.
        D = pdist(Mtry,distfun);

        % If minimal pairwise distance is greater than the current best
        % value, save the coding matrix.
        mindist = min(D);
        if mindist>maxmindist
            M = Mtry;
            maxmindist = mindist;
        end
    end
    
    M(isnan(M)) = 0;
    
    M = clean(M);
    
elseif tf(7) % dense random
    if K<6
        warning(message('stats:designecoc:NotEnoughClassesForDenseRandom'));
    end
    
    % Binary Hamming distance
    distfun = @(y,M) sum( 1 - bsxfun(@times,M,y), 2 )/2;
    
    M = [];
    maxmindist = -Inf;
    
    % Set the number of columns
    [~,L] = log2(K);
    L = 10*L;
    
    % Generate many matrices and choose the one maximizing the minimal
    % Hamming distance between rows. Elements -1 and +1 are generated with
    % probability 1/2 each.
    for n=1:N
        Mtry = rand(K,L);
        iminus = Mtry<1/2;
        iplus = Mtry>1/2;
        Mtry(iminus) = -1;
        Mtry(iplus)  = +1;

        % Compute pairwise distances for all rows and store in vector D.
        D = pdist(Mtry,distfun);
        
        % If minimal pairwise distance is greater than the current best
        % value, save the coding matrix.
        mindist = min(D);
        if mindist>maxmindist
            M = Mtry;
            maxmindist = mindist;
        end
    end
    
    M = clean(M);
    
end

end


function M = clean(M)

L = size(M,2);

badcols = false(1,L);
for i=1:L
    if ~any(M(:,i)==1) || ~any(M(:,i)==-1)
        badcols(i) = true;
        continue;
    end
    
    for j=i+1:L
        if all(M(:,i)==M(:,j)) || all(M(:,i)==-M(:,j))
            badcols(i) = true;
            break;
        end
    end
end

M(:,badcols) = [];
end

