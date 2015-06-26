classdef ChannelStats
    %CHANNELSTATS Summary of this class goes here
    %   Detailed explanation goes here
    
    %drawB( type, threshold )
    %drawT( type, threshold )
    
    %printB
    %printT
    
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
        conditions
        
        tstat
        
        p
        q
    end
    
    methods
        % unique conditions
        function c = get.conditions( obj )
            if ~isempty(obj.variables)
                c = unique(obj.variables.condition); 
            else
                c = [];
            end
        end
        
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
                out = - tinv( abs(pcrit)/2, obj.dfe );
                
            elseif s{1} == 'q'
                t = abs(obj.tstat(:))';
                q = obj.q(:);
                
                out = interp1(q, t, str2num(s{2}));
            end
        end
        
        function out = table( obj )
            beta = obj.beta;
            tstat = obj.tstat;
            p = obj.p;
            q = obj.q;
            se = sqrt(diag(obj.covb));
            dfe = obj.dfe * ones(size(p));
            
            out = [obj.variables table(beta, se, tstat, dfe, p, q)];
        end
        
        stats = ttest( obj, c, b );
        stats = ftest( obj, m );
        stats = jointTest( obj );
        
        draw( obj, vtype, vrange, thresh );
    end
    
    methods (Access = protected)
        newNames = transformNames( obj, T );
    end
  
end

        %% stat calculations
    
%         function A = get.joinTest( obj )
%             [~, ~, idx] = unique( table(obj.variables.source, obj.variables.detector, obj.variables.condition), 'rows', 'stable' );
%             
%             for i = 1:max(idx)
%                 F(i,1) = 
%             end
%             
%         end

