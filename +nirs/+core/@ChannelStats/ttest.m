function [S,haserror] = ttest(obj, c, b, names)
    %% ttest - Tests the null hypothesis c*beta = b
    % 
    % Args:
    %     c - contrast vector/matrix
    %     b - (optional) mean of the null distribution
    %     
    % Example:
    %     % tests the sum and difference of betas
    %     stats.ttest([1 1; 1 -1])


    
    haserror=false;
     if(nargin<3)
            b=[];
            names=[];
     elseif(~isnumeric(b) && nargin<4)
         names=b;
         b=[];
     elseif(nargin<4)
           names=[];
     end
    
     if(length(obj)>1)
         
        for i=1:length(obj)
            [S(i), haserror(i)]=obj(i).ttest(c, b, names);
        end
        S(haserror)=[];
        haserror=any(haserror);
        return
     end
     
     
     if(isstr(c) || iscell(c) || iscellstr(c))
         if(isstr(c)); c=cellstr(c); end;
         if(isstr(names)); names=cellstr(names); end;
         if(isempty(names))
             names=c;
         end
         [c,haserror] = nirs.design.contrastvector(c,obj.conditions,obj.basis);
         
         if(haserror)
            warning(['error processing: ' obj.description]);
            
         end
         
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
    
    nchan = size(obj.probe.link,1);
    
    % sort variables
    [~, icond] = sort(obj.conditions);
    
    
    obj=sorted(obj);
    c = c(:, icond);
    
    % full contrast matrix
    C = kron(eye(nchan), c);

    if nargin < 3 || isempty(b)
        b = zeros(size(C,1),1);
    else
        b = repmat(b(icond), [nchan 1]);
    end
    obj = sorted(obj);
    
    % transform beta
    beta = bsxfun(@minus, C*obj.beta, b);

    % new covariance
    covb = C*obj.covb*C';

    % output
    S = obj;

    S.beta  = beta;
    S.covb  = covb;
    
    if length(S.dfe)>1
        S.dfe = abs(C) * obj.dfe ./ sum(abs(C),2);
    end

    % new condition names
    if nargin < 4 && isempty(names)
        cond = obj.transformNames(c);
    else
        cond = names;
    end
    
    cond = repmat( cond(:), [nchan 1] );

    
        link = repmat( obj.probe.link, [size(c,1) 1] );
        if(isa(S.probe,'nirs.core.ProbeROI'))
            link = nirs.util.sortrows(link, {'ROI', 'type'});
        else
        link = nirs.util.sortrows(link, {'source', 'detector', 'type'});
        end
        S.variables = [link table(cond)];
  
    S.description = 'T-test';
end