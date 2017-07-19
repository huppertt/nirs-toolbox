function [beta,bHat,sigma2]=fitlme(X,y,Z,robust);

if(nargin<4)
    robust=false;
end

N=size(y,1);
q = size(Z,2);

Iq = spdiags(ones(q,1),0,q,q);

w=speye(size(y,1),size(y,1));

D = sqrt(eps(class(X)));
b0 = zeros(size(X,2),1);
b = ones(size(X,2),1);

iter=1;
 while((iter<5) & any(abs(b-b0) > D*max(abs(b),abs(b0))))
    wZ=w*Z;
    wy=w*y;
    wX=w*X;
    
    [R,status,S] = chol(sparse(wZ'*wZ + Iq));
    
    Q1 = ((wX'*wZ)*S) / R;
    R1R1t = wX'*wX - Q1*Q1';
    
    R1 = chol(R1R1t,'lower');
    
    cDeltab = R' \ (S'*((wZ'*wy)));
    %cbeta = R1 \ (X'*y - Q1*cDeltab);
    cbeta = R1 \ (wX'*wy - Q1*cDeltab);
    
    % (9) Compute betaHat and Deltab.
    beta = R1' \ cbeta;
    %Deltab = P1' \ (cDeltab - Q1'*betaHat);
    bHat = S*(R \ (cDeltab - Q1'*beta));
    
    resid=(y - X*beta - Z*bHat);
    
    if(~robust)
        break
    end
    w = wfun(resid);
    b0=b;
    b=beta;
    disp(['Robust fit iteration ' num2str(iter) ' : ' num2str(max(abs(b-b0)))]);
    iter=iter+1;
end


r2 = sum(bHat.^2) + sum(resid.^2);

sigma2 = r2/N;


end


function w = wfun(r)
s = mad(r, 0) / 0.6745;
r = r/s/4.685;

w = (1 - r.^2) .* (r < 1 & r > -1);
w=sparse(diag(w));
end
