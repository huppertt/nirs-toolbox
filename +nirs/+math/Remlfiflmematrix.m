function lme = Remlfiflmematrix(X, beta, Z,PAT);
% Solves the model beta = X*A + Z*B using ReML reguarization

% STore the variables for the end
Z=Z./normest(Z)*normest(X);
%beta=beta-Z*inv(Z'*Z)*Z'*beta;

s=1;
% normalize the model
for i=1:size(Z,2)
    lst=find(Z(:,i)~=0);
    s(i)=max(max(X(lst,:)));
end
s=sqrt(s./norm(s));
for i=1:size(Z,2)
    lst=find(Z(:,i)~=0);
    beta(lst)=beta(lst)/s(i);
    X(lst,:)=X(lst,:)/s(i);
    Z(lst,:)=Z(lst,:)/s(i);
end

s=normest(X);
X=X/s;
beta=beta/s;
Z=Z/s;



% Create the covariance components from PAT and Z
n=size(X,1);
for i=1:size(Z,2)
    Q{i,1}=spdiags(Z(:,i),0,n,n);
end
% Now add the off-diag components
cnt=length(Q);
for i=1:size(PAT,1)
    for j=i+1:size(PAT,2)
        if(PAT(i,j)~=0)
            cnt=cnt+1;
            lst=find(Z(:,i)~=0);
            lst2=find(Z(:,j)~=0);
            
            id=sub2ind([n n],lst,lst2);
            id2=sub2ind([n n],lst2,lst);
            
            s=min(max(Q{i}(:)),max(Q{j}(:)));
            qq=Q{i}+Q{j};
            %qq=(qq~=0)*s;
            qq(id)=PAT(i,j)*s;
            qq(id2)=PAT(i,j)*s;
            Q{cnt,1}=sparse(qq);
        end
    end
end



% NOW SOLVE THE REML MODEL
lambda = ones(length(Q),1);  %Initial guess of lambda. 
tol=1E-4;  %REML goes till tolorence (or max iter) 

tolr=1E-8;  

t=1024;
dF    = Inf;  %Initial decent
cnt=0;  %This is a bookkeeping param for display purposes

maxIter=100;

%% This function was adapted from the SPM function spm_reml.m
%  which is part of the SPM toolkit (http://www.fil.ion.ucl.ac.uk/spm)
%
istruck=0;

for iter=1:maxIter
    Ce = tolr*speye(size(Q{1},1));  %Make sure it stays in numerical precision
    
    % E-step:
    
    %Bound lambda to avoid numerircal prec. issues
    lambda=max(lambda,log(tolr));
    lambda=min(lambda,log(1/tolr));
    
    for i = 1:length(Q)
        Ce = Ce + Q{i}*exp(lambda(i));
    end
    
    %      iCe=blkdiag(inv(Ce(1:end/2,1:end/2)),inv(Ce(1+end/2:end,1+end/2:end)));
    iCe = inv(Ce);
    Xt_iCe = X' * iCe;
    Xt_iCe_X = Xt_iCe * X;
    C_beta_y = inv(Xt_iCe_X);  %Estimate of covariance of beta given the measurements
    
    
    % M-step:
    P = iCe - (iCe*X)*C_beta_y*Xt_iCe;
    PY=P*beta;
    for i=1:size(Q,1)
        PQ_i{i}=P*Q{i};
    end
    
    for i = 1:size(Q,1)
        PQ = PQ_i{i};
        PQt=PQ';
        g(i,1) = -0.5*trace(PQ)*exp(lambda(i)) + norm(0.5*PY'*Q{i}*PY,1)*exp(lambda(i));
        for j = i:size(Q,1)
            PQj = PQ_i{j};
            H(i,j) = -0.5*sum(sum(PQt.*PQj))*exp(lambda(i)+lambda(j));
            H(j,i)=H(i,j);
        end
    end
    
    %Now update the lambda.  dLambda = -inv(H)*g
    I=eye(size(H,1));
    % warning('off','MATLAB:nearlySingularMatrix');
    
    dL = (expm(H*t) - I)*inv(full(H)+eye(size(H))*1E-10)*g;
    %dL= pinv(H)*g;
    df    = g'*dL;
    
    if df > dF - exp(-4) 
        if(istruck<5)
            t = max(2,t/2);    
            istruck+istruck+1;
        else
            break;
        end
    else
        lambda = lambda + dL;
        istruck=0;
    end
   
    dF    = df;
    
    for c=1:cnt, fprintf('\b'); end
    cnt=fprintf('%-5s: %i %5s%e','  ReML Iteration',iter,'...',full(dF));
    
    if dF < tol & iter>5, break; end

    
end

[U,S,V]=svd(full(Ce));
W=(diag(1./sqrt(diag(S)))*(U'+V')/2);

disp('done');
lme=fitlmematrix(W*X, W*beta, [], [], 'FitMethod', 'ML');


return