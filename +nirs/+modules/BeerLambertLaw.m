classdef BeerLambertLaw < nirs.modules.AbstractModule
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        PPF = 5 / 50; % arbitrary scaling factor
        chromophores = {'hbo','hbr'};
    end
    
    methods

        function obj = BeerLambertLaw( prevJob )
           obj.name = 'Beer-Lambert Law';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
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
                    lambda = p.link.type(lst);
                    
                    ext = nirs.media.getspectra( lambda );
                    
                    clist = [1 2]; % hbo and hbr; need to fix this
                    
                    % extinction coefficients
                    E = ext(:,clist);
                    
                    % distances
                    L = p.distances(lst);
                    
                    % mbll model
                    EL = bsxfun( @times, E, L *obj.PPF );
                    iEL = pinv(EL);
                    
                    % calculates chromophore concentration (uM)
                    d(:,lst) = (d(:,lst)*iEL') * 1e6;
                    
                    % new channel type
                    type(lst,1) = obj.chromophores;
                end
                
                p.link.type = type;
                [p.link,idx] = sortrows(p.link,{'type','source','detector'});
                
                data(i).data  = d(:,idx);
                data(i).probe = p;
                
            end
        end
    end
    
end

