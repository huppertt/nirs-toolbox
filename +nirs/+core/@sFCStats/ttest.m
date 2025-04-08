function S = ttest(obj, c, b, names)
    %% ttest - Tests the null hypothesis c*beta = b
    % 
    % Args:
    %     c - contrast vector/matrix
    %     b - (optional) mean of the null distribution
    %     
    % Example:
    %     % tests the sum and difference of betas
    %     stats.ttest([1 1; 1 -1])


    
    
     if(nargin<3)
            b=[];
     end
     if(nargin<4)
           names=[];
     end
    
     if(length(obj)>1)
        for i=1:length(obj)
            S(i)=obj(i).ttest(c, b, names);
        end
         return
     end
     
     
     if(isstr(c) || iscell(c) || iscellstr(c))
         if(isstr(c)); c=cellstr(c); end;
         if(isempty(names))
             names=c;
         end
         c = nirs.design.contrastvector(c,obj.conditions);
         
         %Remove names that didn't exist in this file
         lst=all(c==0,2);
         c(lst,:)=[];
         names={names{~lst}};
     end
     
    if(length(obj)>1)
        for idx=1:length(obj)
            S(idx)=obj(idx).ttest(c,b,names);
        end
        return
    end
    
    c = c';
    if nargin < 3 || isempty(b)
        b = zeros(size(c,2),1);
    end
    b = permute(b,[3 2 1]);
    
    ncon = size(c,2);
    dfe = zeros(1,ncon);
    tbl=obj.table;
    for i=1:length(obj.conditions)
        Z(:,i)=tbl(ismember(tbl.condition,obj.conditions{i}),:).Z;
    end

    if(~isempty(obj.ZstdErr))
        ZstdErr = obj.ZstdErr;
    else
        ZstdErr = [];
    end
    
    for cIdx=1:ncon
        
        Z2(:,cIdx) = Z*c(:,cIdx);
       
        if(~isempty(obj.ZstdErr))
            for i=1:size(ZstdErr,1)
                ZstdErr2(i,cIdx,cIdx) = c(:,cIdx)'*squeeze(ZstdErr(i,:,:))*c(:,cIdx);
            end
        end
        
        dfe(1,cIdx) = sum(bsxfun(@times,abs(c(:,cIdx))',obj.dfe),2) ./ sum(abs(c(:,cIdx)));
        
    end
    
    connections=obj.probe.connections;
    connections=connections(ismember(connections.type,connections.type{1}),:);
    
    cnames=repmat(names,height(connections),1);
    connections=repmat(connections,length(names),1);
    connections.type=cnames;

    S=obj;
    S.probe.connections=connections;
    S.R=tanh(Z2(:));
    S.ZstdErr=ZstdErr2;
    %S.conditions=names;
    S.dfe=dfe;
    
    S.description = 'T-test';
end