function stats = mv_ar_irls( X, Y, Pmax, maxBICSearch)

if nargin < 4 || isempty(maxBICSearch), maxBICSearch = 0; end

if ndims(X) == 3
    [m, n, p] = size(X);
    sameX = false;
else
    [m, p] = size(Y);
    n = size(X,2);
    sameX = true;
    X0 = X;  % keep the original 2D design matrix for ar_irls_fast
    X=repmat(X,[1 1 p]);
end

unstack = @(X) permute( reshape(X, [m n p]), [1 3 2] );
stack	= @(X) blkconst(X); % reshape(permute(X, [1 3 2]),[m p n]));
vec     = @(Y) Y(:);


% initial fit
if sameX
    ss_init = nirs.math.ar_irls_fast(Y, X0, Pmax);
    b = ss_init.beta(:);  % [nCond x nChan] -> column vector
else
    b = [];
    for i = 1:size(Y,2)
        ss_init_i = nirs.math.ar_irls_fast(Y(:,i), X(:,:,i), Pmax);
        b = [b; ss_init_i.beta(:)];
    end
end

b0 = b * 1e16; iter = 0;
mdl=nirs.math.varm(p,2);
EstMdl=[];
nCh = size(Y,2);
found_orders = zeros(nCh, 1);

disp([' Starting model estimate']);
%figure; hold on;
%plot(b); drawnow;
while norm(b-b0)/norm(b0) > .1 && iter < 10
    tic;

    b0 = b;


    Xff=X; r=zeros(size(Y)); Yff=Y;
    bb=reshape(b,[size(X,2),size(Y,2)]);
    prev_orders = found_orders;
    parfor i=1:nCh
            r(:,i)=Y(:,i)-squeeze(X(:,:,i))*bb(:,i);
            if maxBICSearch > 0
                if prev_orders(i) == 0
                    a = nirs.math.ar_fit(r(:,i), min(Pmax, maxBICSearch), false);
                    found_orders(i) = length(a) - 1;
                else
                    a = nirs.math.ar_fit(r(:,i), prev_orders(i), true);
                end
            else
                a = nirs.math.ar_fit(r(:,i), Pmax);
            end
            f = [1; -a(2:end)];
            Xff(:,:,i) = myFilter(f,X(:,:,i));
            Yff(:,i) = myFilter(f,Y(:,i));
    end

    % residual
    % r = vec(Yff) -stack(Xff)*b;
    % r = reshape(r, size(Yff));
    if(~isempty(EstMdl))
        EstMdl = estimate(mdl,r,'Mdl0',EstMdl);
    else
        EstMdl = estimate(mdl,r);
    end
    Yf = infer(EstMdl,Yff);

    for i=1:size(X,2);
        Xf(:,i,:)=infer(EstMdl,squeeze(Xff(:,i,:)));
    end



    %b = mySolve(stack(Xw), vec(Yw));

        
    nCh2 = size(Yf,2);
    b = zeros(n*p, 1);
    r = zeros(size(Yf));
    SS = cell(1, nCh2);
    for i=1:nCh2
        [ss,r(:,i)]=nirs.math.ar_irls_fast(Yf(:,i),Xf(:,:,i),Pmax);
        b((i-1)*n+1:i*n) = ss.beta(:);
        SS{i}=ss;
    end
    %plot(b); drawnow;
    % update iter count
    iter = iter + 1;

    disp(['Iteration-' num2str(iter) ' Time elapsed: ' num2str(toc) 's (delta: ' num2str(norm(b-b0)/norm(b0)) ')']);


end
% 
% b2=reshape(b,[n p]);
% r=[];
% for i=1:size(Y,2);
%     r(:,i)=Yf(:,i)-squeeze(Xf(:,:,i))*b2(:,i); 
% end
% 
% r = vec(Yw) -stack(Xw)*b;
% r = reshape(r, size(Yw));

nCh3 = size(Yf,2);
nCol = size(Xf,2);
w = zeros(size(Yf,1), nCh3);
Xw = zeros(size(Xf));
parfor i=1:nCh3
    wi = SS{i}.w;
    w(:,i) = wi;
    for j=1:nCol
        Xw(:,j,i) = wi .* Xf(:,j,i);
    end
end

stats.beta=reshape(b,[n p]);
stats.sigma2=kron(cov(r),eye(n));
stats.covb = pinv(stack(Xw)'*stack(Xw))*stats.sigma2;

% to make it consistent with the existing GLM, the cov should be
% [cond1[1:nchan],cond2[1:nchan]...].  Currently it is
% [1:ncond]chan1,[1:ncond]chan2...]

covb=zeros(p*n,p*n);
S2=zeros(p*n,p*n);
for i=1:n
    covb([1:p]+(i-1)*p,[1:p]+(i-1)*p)=stats.covb(i:n:end,i:n:end);
    S2([1:p]+(i-1)*p,[1:p]+(i-1)*p)=stats.sigma2(i:n:end,i:n:end);
end
stats.covb=covb;
stats.sigma2=S2;

stats.w=w;
stats.model=EstMdl;

stats.dfe   = min(sum(w,1) - numel(b));

stats.tstat = reshape(b./sqrt(diag(stats.covb)),[n p]);
stats.pval = 2*tcdf(-abs(stats.tstat),stats.dfe);     % two-sided
stats.ppos = tcdf(-stats.tstat,stats.dfe);            % one-sided (positive only)
stats.pneg = tcdf(stats.tstat,stats.dfe);             % one-sided (negative only)

end

function x=blkconst(X)
for i=1:size(X,3)
    x{i}=X(:,:,i);
end
x=blkdiag(x{:});
end

function b = mySolve( X, y )
b = pinv(X'*X)*(X'*y);
end

%%
function out = myFilter( f, y )
    % here we are just making the first value zero before filtering to
    % avoid weird effects introduced by zero padding
    y1 = y(1,:);
    
    y = bsxfun(@minus,y,y1);
    
    out = filter(f, 1, y);
    out = bsxfun(@plus,out,sum(f)*y1); % add the corrected offset back

end

function w = wfun( x )
c = 4.685;
s = mad(x(:), 0) / 0.6745;
w = (1 - (x/s/c).^2) .* (abs(x/s) < c);
end