classdef BeerLambertLaw < nirs.jobs.AbstractJob
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        PPF = 5 / 50; % arbitrary scaling factor that can't realisitically be estimated
        chromophores = {'hbo','hbr'};
    end
    
    methods

        function obj = BeerLambertLaw( prevJob )
           obj.name = 'Beer-Lambert Law';
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
                [~,~,idx] = unique([p.link.source p.link.detector],'rows');
                
                for j = 1:max(idx)
                    lst = idx == j;
                    lambda = p.link.type(lst);
                    
                    ext = nirs.modeling.utilities.getSpectra( lambda );
                    
                    clist = [1 2]; % hbo and hbr; need to fix this
                    
                    % extinction coefficients
                    E = ext(:,clist);
                    
                    % distances
                    L = p.distances(lst);
                    
                    % mbll model
                    EL = bsxfun( @times, E, L *obj.PPF );
                    
                    % calculates chromophore concentration (uM)
                    d(:,lst) = (d(:,lst) / EL) * 1e6;
                    
                    % new channel type
                    type(lst,1) = obj.chromophores;
                end
                
                p.link.type = type;
                
                data(i).data  = d;
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

