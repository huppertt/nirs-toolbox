function stats = HAC_IRLS(Y,X,p,type)

if(nargin<3)
    p=8;
end

if(nargin<4)
    type='HC3';
end

nobs=size(Y,1);
nchan=size(Y,2);
nparam=size(X,2);
beta=zeros(nparam,nchan);

resid=zeros(size(Y));
whac=zeros(size(Y));

tic;
for chanIdx=1:nchan
    
    st=nirs.math.ar_irls(Y(:,chanIdx),X,p);
    fX=filter(st.filter{1},1,X);
    fY=filter(st.filter{1},1,Y(:,chanIdx));
    resid(:,chanIdx)=fY-fX*st.beta;
    
    H = fX*pinv(fX'*fX)*fX';
    h=diag(H);
    whac=computeweights(resid(:,chanIdx),type,h,nobs,nobs-nparam);
    PhiHat=fX'*diag(whac)*fX;
    iXtX=pinv(fX'*fX);
    EstCov=iXtX*PhiHat*iXtX;
    
    % ensure it is pos-definite
    [U,S,V]=svd(EstCov);  
    EstCov=(U*S*U'+V*S*V')/2;
    B{chanIdx} = chol(EstCov);

    beta(:,chanIdx)=st.beta;
    dfe(chanIdx,:)=st.dfe*ones(1,size(X,2)); 
    
end


fprintf( '  Finished.  Time Elapsed %f seconds\n',toc);

Cov = [];
for chanIdx=1:nchan
    C1=[];
    for chanIdx2=1:nchan
        cc = corrcoef(resid(:,[chanIdx chanIdx2]));
        C1 = [C1 B{chanIdx}'*B{chanIdx2}*cc(1,2)];
    end
    Cov=[Cov; C1];
end

Cov=(Cov+Cov')/2;

dfe=reshape(dfe,[],1);

% reorder it so its chan1:end- cond1, chan1:end-cond2...  
% currently its chan1-cond1:end, chan2-cond1:end...
lst=1:size(Cov,1);
lst=reshape(reshape(lst,nparam,nchan)',[],1);
Cov =Cov(lst,lst);


stats.beta = reshape(beta',[],1);

stats.covb=Cov;
stats.se=sqrt(diag(Cov));
stats.tstat=stats.beta./stats.se;
stats.dfe=dfe;
stats.p = 2*tcdf(-abs(stats.tstat),stats.dfe);

return


% -----------------------------
function w = bisqweights(X,r,tune)
% huber bisquare
tiny_s = 1e-6 * std(r);
if tiny_s==0
    tiny_s = 1;
end

[nobs,nparam]=size(X);
[Q,R,perm] = qr(X,0);
if isempty(R) % can only happen if p=0 as we know n>p>=0
    tol = 1;
else
    tol = abs(R(1)) * max(nobs,nparam) * eps(class(R));
end
xrank = sum(abs(diag(R)) > tol);
E = X(:,perm)/R(1:xrank,1:xrank);
h = min(.9999, sum(E.*E,2));
adjfactor = 1 ./ sqrt(1-h);

radj = r .* adjfactor;
s = madsigma(radj,xrank);


wfun = @bisquare;
tune = 4.685;

w = feval(wfun, radj/(max(s,tiny_s)*tune));


return


% --------- weight functions

function w = andrews(r)
r = max(sqrt(eps(class(r))), abs(r));
w = (abs(r)<pi) .* sin(r) ./ r;
return

function w = bisquare(r)
w = (abs(r)<1) .* (1 - r.^2).^2;
return

function w = cauchy(r)
w = 1 ./ (1 + r.^2);
return

function w = fair(r)
w = 1 ./ (1 + abs(r));
return

function w = huber(r)
w = 1 ./ max(1, abs(r));
return

function w = logistic(r)
r = max(sqrt(eps(class(r))), abs(r));
w = tanh(r) ./ r;
return

function w = ols(r)
w = ones(size(r));
return

function w = talwar(r)
w = 1 * (abs(r)<1);
return

function w = welsch(r)
w = exp(-(r.^2));
return


% -----------------------------
function s = madsigma(r,p)
%MADSIGMA    Compute sigma estimate using MAD of residuals from 0
rs = sort(abs(r));
s = median(rs(max(1,p):end)) / 0.6745;
return

% -----------------------------
function w = computeweights(u,type,h,nobs,dfe)

u2 = u.^2;


switch(type)
    case 'CLM'
        
        sse = sum(u2);
        w = repmat(sse/dfe,size(u));
        
    case 'HC0'
        
        w = u2;
        
    case 'HC1'
        
        w = (nobs/dfe)*u2;
        
    case 'HC2'
        
        w = u2./(1-h);
        
    case 'HC3'
        
        w = u2./((1-h).^2);
        
    case 'HC4'
        
        d = min(4,h/mean(h));
        w = u2./((1-h).^d);
    otherwise
        
        error(message('econ:hac:HCWeightsInvalid'))
        
end

return
