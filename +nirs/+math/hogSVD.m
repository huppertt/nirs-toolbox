function [US,V]=hogSVD(L)
% This function computes a higher order generalized SVD
% Do a generalized SVD
% [U1,U2...,V,S1,S2...]=gsvd(L1,L2,...);
% L1 = U1*S1*V'
% L2 = U2*S2*V'

Ls=[];
for i=1:length(L)
    flds=fields(L{i});
    LL=abs(L{i}.(flds{1}))/normest(abs(L{i}.(flds{1})));
    for j=2:length(flds)
        LL=LL+abs(L{i}.(flds{j}))/normest(abs(L{i}.(flds{1})));
        
    end
    Ls=[Ls; LL];
end

[~,~,V]=nirs.math.mysvd(Ls);

%[~,~,V,~,~] = gsvd(Ls(1:31,:),LL);
% A = U*C*X'
% B = V*S*X'
% C'*C + S'*S = I
 

tol=sum(abs(V),2);
V(find(tol<max(tol)*1E-6),:)=0;

US={};
for i=1:length(L)
    US{i}=[];
    flds=fields(L{i});
    for j=1:length(flds)
        US{i}=setfield(US{i},flds{j},L{i}.(flds{j})*V);
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