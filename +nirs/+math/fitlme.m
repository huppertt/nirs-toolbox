function [beta,bHat,sigma2]=fitlme(X,y,Z);

N=size(y,1);
q = size(Z,2);

Iq = spdiags(ones(q,1),0,q,q);

[R,status,S] = chol(Z'*Z + Iq);

Q1 = ((X'*Z)*S) / R;
R1R1t = X'*X - Q1*Q1';

R1 = chol(R1R1t,'lower');

cDeltab = R' \ (S'*((Z'*y)));
%cbeta = R1 \ (X'*y - Q1*cDeltab);
cbeta = R1 \ (X'*y - Q1*cDeltab);

% (9) Compute betaHat and Deltab.
beta = R1' \ cbeta;
%Deltab = P1' \ (cDeltab - Q1'*betaHat);
bHat = S*(R \ (cDeltab - Q1'*beta));


r2 = sum(bHat.^2) + sum((y - X*beta - Z*bHat).^2);

sigma2 = r2/N;
