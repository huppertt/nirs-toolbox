function sttonnica(L,V,mesh,nC)
% Implementation of Spatial-temporal tomograohic NN ICA
% Human Brain Mapping 30:1898?1910 (2009)
% EEG Source Imaging With Spatio-Temporal Tomographic Nonnegative IndependentComponent Analysis
% Pedro A. Valdes-Sosa,Mayrim Vega-Hernandez, Jose Miguel Sanchez-Bornot, Eduardo Mart?nez-Montes, and Mar?´a Antonieta Bobes

lambdasmooth=1;

% Laplacian operator
Lap=nirs.inverse.laplacian(mesh);
LtL=(Lap'*Lap);
LtL=LtL.*(LtL>max(LtL(:))*1E-12);
%Smoothing kernel
I=speye(size(LtL,1),size(LtL,2));
S = inv(I+lambdasmooth*LtL);
S=S.*(S>max(S(:))*1E-12);

[u,s,v]=nirs.math.mysvd(L);
K=u*s;
S=v'*S*v;

% The model is:
% V = K*S*M*G

% Step 1:
% Initialize M with random positive numbers. 
M=rand(size(v,2),nC)*normest(V);

for i=1:100
    
%Step 2:
% Estimate G by ordinary least squares:
X=K*S*M;
G = inv(X'*X)*X'*V;

%Step 3: estimate M
% V = K*S*M*G
X=K*S;
M = inv(X'*X)*X'*V*G'*inv(G*G');

disp(norm(K*S*M*G-V));
end