function [T newTbl] = cp2tform(tbl1,tbl2,rigid)
% This function does an afine tform of all the points in tbl1 to tbl2

if(nargin<3)
    rigid=false;
end

if(isa(tbl1,'table'))
    Pos1=[tbl1.X tbl1.Y tbl1.Z];
    Pos2=[tbl2.X tbl2.Y tbl2.Z];
    [i,j]=ismember(lower(tbl1.Name),lower(tbl2.Name));
else % 3D mats provided
    Pos1=tbl1;
    Pos2=tbl2;
    i=1:size(Pos1,1);
    j=1:size(Pos1,1);
end
Pos1(:,4)=1;
Pos2(:,4)=1;

if(rigid)
    
    

    p=Pos1(find(i),1:3)';
    q=Pos2(j(find(i)),1:3)';
    m = size(p,2);
    n = size(q,2);
    % find data centroid and deviations from centroid
    q_bar = mean(q,2);
    q_mark = q - repmat(q_bar, 1, n);
    
    
    % find data centroid and deviations from centroid
    p_bar = mean(p,2);
    p_mark = p - repmat(p_bar, 1, m);
    
    N = p_mark*transpose(q_mark); % taking points of q in matched order
    
    [U,~,V] = svd(N); % singular value decomposition
    
    R = V*diag([1 1 det(U*V')])*transpose(U);
    
    T = q_bar - R*p_bar;
    
    T=[R T; [0 0 0 1]]';
else
    T = Pos1(find(i),:)\Pos2(j(find(i)),:);
end

if(nargout==2)
    p=Pos1*T;
    newTbl=tbl1;
    newTbl.X=p(:,1);
    newTbl.Y=p(:,2);
    newTbl.Z=p(:,3);
end
    