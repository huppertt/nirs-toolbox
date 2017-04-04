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
    
   
    for i=1:size(obj.Z,1)
        for j=1:size(obj.Z,2)
            Z(i,j,:)=c*squeeze(obj.Z(i,j,:));
        end
    end
    
    S=obj;
    S.R=tanh(Z);
    S.conditions=names;
    S.dfe=mean(bsxfun(@times,abs(c),obj.dfe));
    
    S.description = 'T-test';
end