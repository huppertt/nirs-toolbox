classdef BeerLambertLaw < nirs.modules.AbstractModule
%% BeerLambertLaw - Converts optical density to hemoglobin.
% 
% dOD(lambda) = ext(hbo, lambda) * conc(hbo) * distance * PPF + ...
%         ext(hbr, lambda) * conc(hbr) * distance * PPF;
%
% Options: 
%     tune - number of standard deviations to define an outlier
    
    properties
        PPF = 5 / 50;   % partial pathlength factor 
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
                    
                    assert( length(lst) > 1 )
                    
                    lambda = p.link.type(lst);
                    
                    ext = nirs.media.getspectra( lambda );
                    
                    clist = [1 2]; % hbo and hbr; need to fix this
                    
                    % extinction coefficients
                    E = ext(:,clist);
                    
                    % distances
                    L = p.distances(lst);
                    
                    % mbll model
                    EL = bsxfun( @times, E, L*obj.PPF );
                    iEL = pinv(EL);
                    
                    % calculates chromophore concentration (uM)
                    d(:,lst) = (d(:,lst)*iEL') * 1e6;
                    
                    % new channel type
                    type(lst,1) = {'hbo', 'hbr'};
                end
                
                p.link.type = type;
                [p.link,idx] = sortrows(p.link,{'source','detector','type'});
                
                data(i).data  = d(:,idx);
                data(i).probe = p;
                
            end
        end
    end
    
end

