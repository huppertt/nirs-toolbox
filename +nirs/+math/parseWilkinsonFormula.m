function [X, Z, names] = parseWilkinsonFormula( formula, tbl, iscentered, dummycoding )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    if nargin < 3
        iscentered = true;
    end
    
    if nargin < 4
        dummycoding = 'reference';
    end

    X = []; Z = []; names = {};
    
    %% split dependent and independent
    terms = strtrim( strsplit(formula,'~') );
    
    lhs = terms{1};
    rhs = terms{2};

    % check for unsupported operators
    if length( strsplit(rhs,{':','^'}) ) > 1
        error('This only supports ''*'' operator currently.')
    end
    
    % split additive terms
    rhs = strtrim( strsplit(rhs,'+') );

    %% constant term
    lst = strcmp(rhs,'-1');
    
    if any(lst)
        isconstant = false;
    else
        isconstant = true;
    end
    rhs = rhs(~lst);
    
    if isconstant
        lst = strcmp(rhs,'1');
        rhs = rhs(~lst);
        X = [X; ones(size(tbl,1),1)];
        names{end+1,1} = 'intercept';
    end
    
    %% additional terms
    for i = 1:length(rhs)
        %% RFX terms
        if strcmp( rhs{i}([1 end]), '()' )
            T = strtrim( strsplit( rhs{i}(2:end-1),'|' ) );
            if ~strcmp(T{1},'1')
                error( 'This only supports random intercept models.' )
            end
            
            T = strtrim( strsplit( T{2},'*' ) );
            
            isnum = true;
            iscat = true;
            for j = 1:length(T)
                isnum = isnum & isnumeric( tbl.(T{j}) );
                iscat = iscat & iscellstr( tbl.(T{j}) );
            end
            
            assert( ((isnum && iscat) == false) && iscat )
            
            z = parseCategorical( T, tbl );
            Z = [Z z];
            
        
        %% FFX term    
        else 
            T = strtrim( strsplit( rhs{i},'*' )  );
            
%             x = []; n = {};
%             
%             % parse each piece
%             for j = 1:length( T )
%                 
%                 isnum = isnumeric( tbl.(T{j}) );
%                 iscat = iscellstr( tbl.(T{j}) );
%                 
%                 assert( isnum || iscat )
%                 
%                 if iscat
%                     [a, b] = parseCategorical( T{j}, tbl );
%                     if strcmp(dummycoding,'reference')
%                         x = [x a(:,2:end)];
%                         n = [n b(2:end)];
%                     elseif strcmp(dummycoding, 'full' )
%                         x = [x a];
%                         n = [n; b];
%                     end
%                 else
%                     if iscentered
%                         x = [x (tbl.(T{j}) - mean(tbl.(T{j})))];
%                     else
%                         x = [x tbl.(T{j})];
%                     end
%                 	n = [n; T{j}];
%                 end
%             end
%             
%             % pair-wise multiplication
%             for j = 1
%                 
%             end
%             
%             %%% WRITE PAIRWISE MULTIPLY AND CALL RECURSIVELY
            
            isnum = true;
            iscat = true;
            for j = 1:length(T)
                isnum = isnum & isnumeric( tbl.(T{j}) );
                iscat = iscat & iscellstr( tbl.(T{j}) );
            end
            
            assert( (isnum & iscat) == false )
            
            if iscat
                [x,n] = parseCategorical( T, tbl );
                
                if strcmp(dummycoding,'reference')
                    x = x(:,2:end);
                    n = n(2:end);
                end

                X = [X x];
                names = [names; n];
            else
                x = ones(size(tbl,1),1);
                for j = 1:length(T)
                    if iscentered
                        x = x .* (tbl.(T{j}) - mean(tbl.(T{j})));
                    else
                        x = x .* tbl.(T{j});
                    end
                end
                
                X = [X x];
                names = [names; strjoin(T,':')];
            end
        end
    end
    
%     %% Dependent Variable
%     y = tbl.(lhs);

end

function [X, names] = parseCategorical( T, tbl )
        if ischar(T)
            T = {T};
        end

        n = {};
        for j = 1:length(T)
            n = [n tbl.(T{j})];
        end

        names = {};
        for j = 1:size(n,1)
            names{j,1} = strjoin(n(j,:),':');
        end

        [names,~,I] = unique(names,'stable'); 
        
        % if you input strings to dummyvar it will sort them
        X = dummyvar( I );
end

