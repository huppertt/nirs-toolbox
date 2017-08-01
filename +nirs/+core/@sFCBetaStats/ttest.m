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
    beta = zeros([size(obj.beta,1) size(obj.beta,2) ncon]);
    covb = zeros([size(obj.beta,1) size(obj.beta,2) ncon ncon]);
    
    for cIdx = 1:size(c,2)
        
        dfe(1,cIdx) = sum(bsxfun(@times,abs(c(:,cIdx))',obj.dfe),2) ./ sum(abs(c(:,cIdx)));
        
        tmpC = permute(c(:,cIdx),[3 2 1]);
        
        beta(:,:,cIdx) = sum( bsxfun( @times , obj.beta , tmpC ) , 3);
        
        % new covariance
        tmpCT = permute(tmpC,[1 2 4 3]);
        covb(:,:,cIdx,cIdx) = sum(sum( bsxfun( @times , bsxfun( @times , obj.covb , tmpC ) , tmpCT ) ,4),3);

    end
    
    beta = bsxfun( @minus , beta , b );
    
    % output
    S = nirs.core.sFCStats;
    S.type = obj.type;
    S.probe = obj.probe;
    S.demographics = obj.demographics;
    S.conditions=names;
    S.R = tanh(beta);
    S.ZstdErr = sqrt(covb);
    S.dfe = dfe;   
    S.description = 'T-test';
end