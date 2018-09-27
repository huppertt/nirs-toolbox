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
                if(isa(data,'nirs.core.Data'))
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
              lstrm=[];
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
                    lst2=find(lst);
                    if(length(lst2)>2)
                        lstrm=[lstrm lst2(3:end)];
                        lst2=lst2(1:2);
                    end
                    d(:,lst2) = (d(:,lst)*iEL') * 1e6;
                    
                    % new channel type
                    type(lst2,1) = {'hbo', 'hbr'};
                end
                type(lstrm(find(lstrm<=length(type))))=[];
                d(:,lstrm)=[];
                p.link(lstrm,:)=[];
                p.link.type = type;
                
                p.link.type = type;
                
                if(~ismember('source',p.link.Properties.VariableNames) & ...
                        ismember('ROI',p.link.Properties.VariableNames))
                    [p.link,idx] = nirs.util.sortrows(p.link,{'ROI','type'});
                else
                    [p.link,idx] = nirs.util.sortrows(p.link,{'source','detector','type'});
                end
                
                data(i).data  = d(:,idx);
                data(i).probe = p;
                elseif(isa(data,'nirs.core.ChannelStats'))
                    
                    
                    
                    if(~ismember('source',data(i).probe.link.Properties.VariableNames) & ...
                            ismember('ROI',data(i).probe.link.Properties.VariableNames))
                       data(i)=data(i).sorted({'ROI','type'});
                    else
                        % sort channels
                        data(i)=data(i).sorted({'source','detector','type'});
                    end
                    var=data(i).variables;
                    var.type=[];
                    [~,idx]=unique(var);
                    
                    var2=data(i).variables;
                    var2.cond=[];
                    var2=unique(var2);
                    
                    clear type;
                    
                    F=zeros(height(data(i).variables));
                    N=cell(height(data(i).variables),1);
                    N2=cell(height(data(i).probe.link),1);
                    for j = 1:length(idx)
                        lst = find(ismember(var,var(idx(j),:)));
                        lst2 = find(ismember(var2(:,1:2),var(idx(j),1:2)));
                        assert( length(lst) > 1 )
                        
                        lambda = data(i).variables.type(lst);
                        ext = nirs.media.getspectra( lambda );
                        
                        clist = [1 2]; % hbo and hbr; need to fix this
                        
                        % extinction coefficients
                        E = ext(:,clist);
                        
                        % distances
                        L = data(i).probe.distances(lst2);
                        L=max(L,1);  % avoid issues with the short (0) seperation values
                        
                        % mbll model
                        EL = bsxfun( @times, E, L*obj.PPF );
                        iEL = pinv(EL);
                        
                        % calculates chromophore concentration (uM)
                        F(lst,lst)=iEL' * 1e6;
                        N(lst)={'hbo', 'hbr'};
                        N2(lst2)={'hbo', 'hbr'};
                        % new channel type
                        
                    end
                
                    data(i).beta=F*data(i).beta;
                    data(i).covb=F*data(i).covb*F';
                    data(i).variables.type=N;
                    data(i).probe.link.type=N2;
                    
                    
                else
                    error('unsupported type');
                end
            end
        end
    end
    
end

