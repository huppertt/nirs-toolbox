function [US,V]=hogSVD(L)
% This function computes a higher order generalized SVD
% Do a generalized SVD
% [U1,U2...,V,S1,S2...]=gsvd(L1,L2,...);
% L1 = U1*S1*V'
% L2 = U2*S2*V'

fld={};
for i=1:length(L)
    f=fields(L{i});
    for j=1:length(f)
        fld{end+1,1}=f{j};
    end
end


L2=[]; cnt=1;
for i=1:length(L)
   f=fields(L{i});
   L3=[];
   for j=1:length(fld)
         id=find(ismember(f,fld{j}));
         if(~isempty(id))
            L3=[L3 L{i}.(f{id})];
         else
             L3=[L3 0*L{i}.(f{1})];
         end
    end 
    
   L2=[L2; L3];
end

% % I am using parafac instead of SVD to allow for missing (NaN) data
% F=parafac(L2,size(L2,1),[1e-6 1 0 1 NaN 50],1);
% V=F{2};
% V=orth(V);

if(length(fld)>1)
    LL=reshape(L2,size(L2,1),[],length(fld));
    for i=1:length(fld)
        LL(:,:,i)=LL(:,:,i)/norm(LL(:,:,i));
    end
    a=squeeze(sum(sum(abs(LL),3),1));
    lst=find(a>max(a)/1000);
    
    a=parafac(LL(:,lst,:),size(LL,1));
    V=zeros(size(LL,2),size(a{2},2));
    V(lst,:)=a{2};
else
    [U,S,V]=nirs.math.mysvd(L2);
end

% tol=sum(abs(V),2);
% V(find(tol<max(tol)*1E-6),:)=0;

US={};
for i=1:length(L)
    US{i}=struct;
    f=fields(L{i});
    L3=[];
    for j=1:length(fld)
         id=find(ismember(f,fld{j}));
         if(~isempty(id))
            US{i}=setfield(US{i},f{id}, L{i}.(f{id})*V);
         else
             US{i}=setfield(US{i},fld{j}, 0*L{i}.(f{1})*V);
         end
     
    end 
     
end


return

% function [US,V]=hogSVD(L)
% % This function computes a higher order generalized SVD
% % Do a generalized SVD
% % [U1,U2...,V,S1,S2...]=gsvd(L1,L2,...);
% % L1 = U1*S1*V'
% % L2 = U2*S2*V'
% 
% 
% L2=[]; cnt=1;
% for i=1:length(L)
%    f=fields(L{i});
%    L3=[];
%     for j=1:length(f)
%      %   L2(:,:,cnt)=sparse(L{i}.(f{j}));
%      %   cnt=cnt+1;
%          L3=[L3; L{i}.(f{j})];
%          m(i)=size(L{i}.(f{j}),1);
%     end 
%    L2=[L2; L3];
% end
% 
% % % I am using parafac instead of SVD to allow for missing (NaN) data
% % F=parafac(L2,size(L2,1),[1e-6 1 0 1 NaN 50],1);
% % V=F{2};
% % V=orth(V);
% 
% 
% [U,S,V]=nirs.math.mysvd(L2);
% 
% 
% % tol=sum(abs(V),2);
% % V(find(tol<max(tol)*1E-6),:)=0;
% 
% US={};
% for i=1:length(L)
%     US{i}=struct;
%     f=fields(L{i});
%     for j=1:length(f)
%         US{i}=setfield(US{i},f{j},L{i}.(f{j})*V);   
%     end
%     
% end
% 
% 
% 
% return