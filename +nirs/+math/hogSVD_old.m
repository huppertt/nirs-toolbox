function [US,V]=hogSVD(L)
% This function computes a higher order generalized SVD
% Do a generalized SVD
% [U1,U2...,V,S1,S2...]=gsvd(L1,L2,...);
% L1 = U1*S1*V'
% L2 = U2*S2*V'

Ls=abs(L{1,1});
for i=2:size(L,1)
    for j=1:size(L,2)
        Ls=Ls+abs(L{i,j});
    end
end

[~,~,V]=nirs.math.mysvd(Ls);


tol=sum(abs(V),2);
V(find(tol<max(tol)*1E-6),:)=0;

for i=1:size(L,1)
    for j=1:size(L,2)
        US{i,j}=L{i,j}*V;
    end
end


% 
% %deal with the trivial cases first
% if(length(L)==1)
%     [U{1},S{1},V]=nirs.math.mysvd(full(L{1}));
%     return;
% elseif(length(L)==2)
%     [U{1},U{2},V,S{1},S{2}] = gsvd(full(L{1}),full(L{2}),0);
%     return;
% else
%     error('not implmented yet');
% end

%Now the harder case of multiple pairs


return