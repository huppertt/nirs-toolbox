classdef ChannelFStats
    %ANOVASTATS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        variables	% variable names
        F           % F-stats
        df1         % degrees of freedom 1
        df2         % degrees of freedom 2
        probe       % probe geometry
    end
    
    properties ( Dependent = true )
        p
        q
    end
    
    methods
        % critical value
        function fcrit = getCritF( obj )
            fcrit = zeros(length(obj.names),1);
            for i = 1:length(obj.names)
                fcrit(i) = 1/finv( obj.pcrit, obj.df2(i), obj.df1(i) );
            end
        end
        
        % p values
        function p = get.p( obj )
            p = zeros(size(obj.F));
            for i = 1:numel(obj.F)
                p(i) = fcdf( 1./obj.F(i), obj.df2(i), obj.df1(i) );
            end
        end
        
        % q values
        function q = get.q( obj )
            q = reshape( nirs.math.fdr( obj.p(:) )', size(obj.p) );
        end
        
        % stats table
        function out = table( obj )
            F = obj.F;
            p = obj.p;
            q = obj.q;
            df1 = obj.df1;
            df2 = obj.df2;
            
            out = [obj.variables table(F, df1, df2, p, q)];
        end
        
        % draw
        draw( obj,  frange, idx );
    end
end