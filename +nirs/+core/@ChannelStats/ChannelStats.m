classdef ChannelStats
    %CHANNELSTATS Summary of this class goes here
    %   Detailed explanation goes here
    
    %drawB( type, threshold )
    %drawF( threshold )
    %drawT( type, threshold )
    
    %printB
    %printF
    %printT
    
    %F -- dependent
    
    properties
        description
        
        variables
        
        beta            % nconds x nchannels
        covb            % nconds x nconds x nchannels
        dfe          	% degrees of freedom
        
        probe           % probe geometry
        demographics    % subject demographics
    end
    
    properties ( Dependent = true )
        tstat
        
        p
        q
    end
    
    methods
        % t statistic calculation
        function tstat = get.tstat( obj )
            tstat = obj.beta ./ sqrt(diag(obj.covb));
        end
        
        % p value calculation
        function p = get.p( obj )
            t = obj.tstat;
            p = 2*tcdf(-abs(t), obj.dfe);
        end
        
        % q values
        function q = get.q( obj )
            q = reshape( nirs.math.fdr( obj.p(:) )', size(obj.p) );
        end
        
        % critical value
        function out = getCritT( obj, s )
            % takes in string in form of 'p < 0.05' or 'q < 0.10'
            s = strtrim( strsplit( s, '<' ) );
            
            if s{1} == 'p'
                pcrit = str2num(s{2});
                out = - tinv( abs( pcrit )/2, obj.dfe );
                
            elseif s{1} == 'q'
                t = abs(obj.tstat(:))';
                q = obj.q(:);
                
                out = interp1(q, t, str2num(s{2}));
            end
        end
    
        % linear transform of variables
        function S = ttest(obj, C, b)
            % test hypothesis that C*beta = b
            
            if nargin < 3
                b = zeros(size(C,1),1);
            end
            
            % transform beta
            beta = bsxfun(@minus, C*obj.beta, b);
            
            % calculate new covariance
            for i = 1:size(beta,2)
                covb(:,:,i) = C*obj.covb(:,:,i)*C';
                tstat(:,i)  = beta(:,i) ./ sqrt( diag(covb(:,:,i)) );
            end
            
            S = obj;
            
            S.beta  = beta;
            S.covb  = covb;
            S.names = obj.transformNames(C);
        end
        
        function S = ftest(obj, m)
            % this uses Hotelling's T squared test for joint 
            % hypothesis testing
            
            if nargin == 1
                m = ones(size(obj.beta,1),1) > 0;
            end
            
            if ~islogical(m)
                m = m > 0;
                warning('Converting mask to true/false.')
            end
                        
            n = obj.dfe;
            k = size(obj.beta,1);
            
            for i = 1:size(m,1)
                for j = 1:size(obj.beta,2)
                    n = obj.dfe;
                    k = sum(m(i,:));
                    
                    b = obj.beta(m(i,:),j);
                    T2(i,j) = b'*pinv(obj.covb(m(i,:),m(i,:),j))*b;
                    
                    F(i,j)      = (n-k) / k / (n-1) * T2(i,j);
                    df1(i,j)    = k;
                    df2(i,j)    = n-k;
                end
                names{i,1} = ['F' num2str(i)];
            end
            
            S = nirs.core.AnovaStats();
            
            S.names = names;
            S.F   = F;
            S.df1 = df1;
            S.df2 = df2;
            S.probe = obj.probe;
        end
        
        function draw( obj, vtype, vrange, vcrit )
            % type is either beta or tstat
            if nargin < 2, vtype = 'tstat'; end
            
            values = obj.(vtype);
            
            % range to show
            if nargin < 3
                vmin = min( values(:) );
                vmax = max( values(:) );
                vrange = [vmin vmax];
            end
            
            if nargin < 4, vcrit = 0; end
            
            % loop through var names
            h = []; % handles
            types = obj.probe.link.type;
            
            if any(isnumeric(types))
                types = cellfun(@(x) {num2str(x)}, num2cell(types));
            end
            
            utypes = unique(types, 'stable');
            
            for iName = 1:length( obj.names )
                
                for iType = 1%:length(utypes)
                    lst = strcmp( types, utypes(iType) );
                    
                    h(end+1) = figure;
                    v = values(iName, lst);
                    obj.probe.draw( v, vrange, vcrit );
                    title([utypes(iType) ' : ' obj.names{iName}], 'Interpreter','none')
                end
            end
        end
    end
    
    methods (Access = protected)
        % this function generates new names based on a linear tranformation
        % T of the variables
        function newNames = transformNames( obj, T )
            names = obj.names;
            for i = 1:size(T,1)
                newNames{i} = '';
                for j = 1:size(T,2)
                    c = T(i,j);
                    if c == 1
                        newNames{i} = [newNames{i} '+' names{j}];
                    elseif c == -1
                        newNames{i} = [newNames{i} '-' names{j}];
                    elseif c > 0
                        newNames{i} = [newNames{i} '+' num2str(c) names{j}];
                    elseif c < 0
                        newNames{i} = [newNames{i} num2str(c) names{j}];
                    end
                end
                
                if newNames{i}(1) == '+'
                    newNames{i} = newNames{i}(2:end);
                end
            end
            
            newNames = newNames(:);
        end
    end
  
end

