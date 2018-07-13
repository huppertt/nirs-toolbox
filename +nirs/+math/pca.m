function [PCs,projmat,lambda] = pca( X , nPC )
if nargin<2, nPC = size(X,2); end

% Remove mean
X = bsxfun(@minus,X,mean(X));

% Compute covariance
C = nancov(X);

% Compute eigenvectors and eigenvalues
[V,D] = eig(C,'vector');
[D,sort_ind] = sort(D,'descend');
V = V(:,sort_ind);

% Find smallest eigenvalue needed
if 0<nPC && nPC<1
    g = cumsum(D) / sum(D);
    lastInd = find(g>=nPC,1,'first');
elseif 1<=nPC && nPC<=length(D) && nPC==round(nPC)
    lastInd = nPC;
else
    error('Bad number of PCs: %g',nPC);
end

% Reduce dimension
V = V(:,1:lastInd);
D = D(1:lastInd);

% Return outputs
PCs = X * V;
projmat = V;
lambda = D;

end