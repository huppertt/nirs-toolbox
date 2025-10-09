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
    Z = zeros([size(obj.R,1) size(obj.R,2) ncon]);
    if(~isempty(obj.ZstdErr))
        ZstdErr = zeros([size(obj.R,1) size(obj.R,2) ncon ncon]);
    else
        ZstdErr = [];
    end
    
    for cIdx=1:ncon
        
        tmpC = permute(c(:,cIdx),[3 2 1]);
        
        Z(:,:,cIdx) = sum( bsxfun( @times , obj.Z , tmpC ) , 3);
       
        if(~isempty(obj.ZstdErr))

            tmpCT = permute(tmpC,[1 2 4 3]);
            ZstdErr(:,:,cIdx,cIdx) = sum(sum( bsxfun( @times , bsxfun( @times , obj.ZstdErr , tmpC ) , tmpCT ) ,4),3);

        end
        
        dfe(1,cIdx) = sum(bsxfun(@times,abs(c(:,cIdx))',obj.dfe),2) ./ sum(abs(c(:,cIdx)));
        
    end
    
    Z = bsxfun( @minus , Z , b );
    
    S=obj;
    S.R=tanh(Z);
    S.ZstdErr=ZstdErr;
    S.conditions=names;
    S.dfe=dfe;
    
    S.description = 'T-test';

    if(~isempty(obj.pvalue_fixed))
        if(all(ismember(names,obj.conditions)))
            % the conditions are directly a subset of the originals and we
            % can keep the fixed pvalues
            mask=zeros(size(obj.t));
            mask(:,:,ismember(obj.conditions,names))=1;
            if(isa(obj.pvalue_fixed,'nirs.bootstrapping.bootstrap_result'))
                lst=find(~mask);
                S.pvalue_fixed.truth(lst)=[];
                S.pvalue_fixed.ecdf(:,lst)=[];
                S.pvalue_fixed.value_bins(:,lst)=[];
            else
                S.pvalue_fixed=obj.pvalue_fixed(find(mask));
            end
        else
            warning('New conditions invalidates permutation p-values');
            S.pvalue_fixed=[];
        end


end