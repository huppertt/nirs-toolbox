classdef ReverseBeerLambert < nirs.functional.AbstractModule
% This converts from concentration to optical density.  This is not really
% useful for anything except simulation.
    
    properties
        PPF = 5 / 50; % arbitrary scaling factor
        lambda = [690 830]';
    end
    
    methods

        function obj = ReverseBeerLambert( prevJob )
           obj.name = 'Convert Hb to OD';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = execute( obj, data )
            for i = 1:length(data)
                d = data(i).data;
                p = data(i).probe;
                
                % sort channels
                [p.link, idx] = sortrows(p.link,{'source','detector','type'});
                d = d(:,idx);
                
                % unique source-detector pairs
                [~,~,idx] = unique([p.link.source p.link.detector],'rows','stable');
                
                for j = 1:max(idx)
                    lst = idx == j;
                    
                    ext = nirs.getSpectra( obj.lambda );
                    
                    clist = [1 2]; % hbo and hbr; need to fix this
                    
                    % extinction coefficients
                    E = ext(:,clist);
                    
                    % distances
                    L = p.distances(lst);
                    
                    % mbll model
                    EL = bsxfun( @times, E, L *obj.PPF );
                    
                    % calculates chromophore concentration (uM)
                    d(:,lst) = (d(:,lst)*EL') * 1e-6;
                    
                    % new channel type
                    type(lst,1) = obj.lambda;
                end
                
                p.link.type = type;
                [p.link,idx] = sortrows(p.link,{'type','source','detector'});
                
                data(i).data  = d(:,idx);
                data(i).probe = p;
                
            end
        end
        
        function options = getOptions( obj )
            options = [];
        end
           
        function obj = putOptions( obj, options )
        end
        
    end
    
end

