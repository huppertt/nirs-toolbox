classdef AnovaStats
    %ANOVASTATS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        names	% variable names
        F       % F-stats
        df1     % degrees of freedom 1
        df2     % degrees of freedom 2
        probe   % probe geometry
        
        pcrit = 0.05;
    end
    
    properties ( Dependent = true )
        Fcrit
        p
    end
    
    methods
        function fcrit = get.Fcrit( obj )
            fcrit = zeros(length(obj.names),1);
            for i = 1:length(obj.names)
                fcrit(i) = finv( 1-obj.pcrit, obj.df1(i), obj.df2(i) );
            end
        end
        
        function p = get.p( obj )
            p = zeros(size(obj.F));
            for i = 1:length(obj.names)
                p(i,:) = fcdf( 1./obj.F(i,:), obj.df2(i), obj.df1(i) );
            end
        end
        
        function draw( obj, frange, idx )
            
            F       = obj.F;
            fcrit   = obj.Fcrit;
            
            if nargin < 2
                frange = [0 ceil(max(abs(F(:))))];
            end
                        
            % loop through var names
            h = []; % handles
            types = obj.probe.link.type;
            
            if any(isnumeric(types))
                types = cellfun(@(x){num2str(x)}, num2cell(types));
            end
            
            utypes = unique(types, 'stable');
            
            for iName = 1:length( obj.names )
                
                for iType = 1:length(utypes)
                    lst = strcmp( types, utypes(iType) );
                    
                    h(end+1) = figure;
                    f = F(iName, lst);
                    obj.probe.draw( f, frange, fcrit(iName) );
                    title([utypes(iType) ' : ' obj.names{iName}], 'Interpreter','none')
                end
            end
        end
    end
    
end