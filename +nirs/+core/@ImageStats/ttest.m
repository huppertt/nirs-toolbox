function S = ttest(obj, c, b,names)
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
    
    nchan = length(obj.beta)/length(obj.conditions);

    % sort variables
    [~, icond] = sort(obj.conditions);
    obj = obj.sorted();
    
    
    
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
    
    
    c = c(:, icond);
    
    % full contrast matrix
    C = kron(speye(nchan,nchan), c);

    if nargin < 3
        b = zeros(size(C,1),1);
    else
        b = repmat(b(icond), [nchan 1]);
    end

    % transform beta
    beta = bsxfun(@minus, C*obj.beta, b);

    % new covariance
    covb=C*obj.covb_chol;
  
    % output
    S = obj;

    S.beta  = beta;
    S.covb_chol  = covb;
    S.typeII_StdE = C*obj.typeII_StdE;
    
    % new condition names
    if(isempty(names))
        cond = obj.transformNames(c);
    else
        cond=names;
    end
    cond = repmat( cond(:), [nchan 1] );

    var = obj.variables;
    type=unique(var.cond);
    lst=find(ismember(var.cond,type{1}));
    var=var(lst,:);
    var=repmat(var,size(c,1),1);
    var.cond=cond;
    S.variables = var;
end