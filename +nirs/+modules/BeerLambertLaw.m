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
           obj.citation='Extinction coef from: Jacques, Steven L. "Optical properties of biological tissues: a review." Physics in medicine and biology 58.11 (2013): R37.';
           
        end
        
        function data = runThis( obj, data )
            for i = 1:numel(data)
                d = data(i).data;
                p = data(i).probe;
                
                
                if(~ismember('source',p.link.Properties.VariableNames) & ...
                        ismember('ROI',p.link.Properties.VariableNames))
                    [p.link, idx] = nirs.util.sortrows(p.link,{'ROI','type'});
                    d = d(:,idx);
                    
                    % unique source-detector pairs
                    [~,~,idx] = nirs.util.uniquerows(table([p.link.ROI]));
                else
                    % sort channels
                    [p.link, idx] = nirs.util.sortrows(p.link,{'source','detector','type'});
                    d = d(:,idx);
                    
                    % unique source-detector pairs
                    [~,~,idx] = nirs.util.uniquerows(table([p.link.source p.link.detector]));
                end
                
                
               clear type;
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
                    L=max(L,1);  % avoid issues with the short (0) seperation values
                    
                    % mbll model
                    EL = bsxfun( @times, E, L*obj.PPF );
                    iEL = pinv(EL);
                    
                    % calculates chromophore concentration (uM)
                    d(:,lst) = (d(:,lst)*iEL') * 1e6;
                    
                    % new channel type
                    type(lst,1) = {'hbo', 'hbr'};
                end
                
                p.link.type = type;
                
                if(~ismember('source',p.link.Properties.VariableNames) & ...
                        ismember('ROI',p.link.Properties.VariableNames))
                    [p.link,idx] = nirs.util.sortrows(p.link,{'ROI','type'});
                else
                    [p.link,idx] = nirs.util.sortrows(p.link,{'source','detector','type'});
                end
                
                data(i).data  = d(:,idx);
                data(i).probe = p;
                
            end
        end
    end
    
end

