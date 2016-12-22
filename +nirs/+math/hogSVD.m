function [US,V]=hogSVD(L)
% This function computes a higher order generalized SVD
% Do a generalized SVD
% [U1,U2...,V,S1,S2...]=gsvd(L1,L2,...);
% L1 = U1*S1*V'
% L2 = U2*S2*V'


L2=[];
for i=1:length(L)
    f=fields(L{i});
   for j=1:length(f)
        L2=[L2; L{i}.(f{j})];
        m(i)=size(L{i}.(f{j}),1);
   end 
end

[U,S,V]=nirs.math.mysvd(L2);
% % I am using parafac instead of SVD to allow for missing (NaN) data
% Factor=parafac(L2,sum(m),[1 1]);
% V=Factor{2};

% tol=sum(abs(V),2);
% V(find(tol<max(tol)*1E-6),:)=0;

US={};
for i=1:length(L)
    US{i}=[];
    flds=fields(L{i});
    for j=1:length(flds)
        lst=find(ismember(flds,f{j}));
        US{i}=setfield(US{i},flds{j},L{i}.(flds{j})*V);
    end
end



return