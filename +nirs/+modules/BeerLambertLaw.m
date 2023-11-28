classdef BeerLambertLaw < nirs.modules.AbstractModule
%% BeerLambertLaw - Converts optical density to hemoglobin.
% 
% dOD(lambda) = ext(hbo, lambda) * conc(hbo) * distance * PPF + ...
%         ext(hbr, lambda) * conc(hbr) * distance * PPF;
%
    
    properties
        PPF = 5 / 50;   % partial pathlength factor 
       % PPF = @(lambda,data)nirs.media.frontal_DPF_model(lambda,data,'age');
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
                    if(~isempty(p.fixeddistances))
                        p.fixeddistances=p.fixeddistances(idx);
                    end
                    % unique source-detector pairs
                    [~,~,idx] = nirs.util.uniquerows(table([p.link.ROI]));
                else
                    % sort channels
                    [p.link, idx] = nirs.util.sortrows(p.link,{'source','detector','type'});
                    d = d(:,idx);
                    if(~isempty(p.fixeddistances))
                        p.fixeddistances=p.fixeddistances(idx);
                    end
                    % unique source-detector pairs
                    [~,~,idx] = nirs.util.uniquerows(table([p.link.source p.link.detector]));
                end
                
                
               clear type;
%               lstrm=[];
%                 for j = 1:max(idx)
%                     lst = idx == j;
%                     
%                     assert( length(lst) > 1 )
%                     
%                     lambda = p.link.type(lst);
%                     ext = nirs.media.getspectra( lambda );
%                     
%                     clist = [1 2]; % hbo and hbr; need to fix this
%                     
%                     % extinction coefficients
%                     E = ext(:,clist);
%                     
%                     % distances
%                     L = p.distances(lst);
%                     L=max(L,1);  % avoid issues with the short (0) seperation values
%                     
%                     % mbll model
%                     EL = bsxfun( @times, E, L*obj.PPF );
%                     iEL = pinv(EL);
%                     
%                     % calculates chromophore concentration (uM)
%                     lst2=find(lst);
%                     if(length(lst2)>2)
%                         lstrm=[lstrm lst2(3:end)];
%                         lst2=lst2(1:2);
%                     end
%                     d(:,lst2) = (d(:,lst)*iEL') * 1e6;
%                     
%                     % new channel type
%                     type(lst2,1) = {'hbo', 'hbr'};
%                 end
%                 type(lstrm(find(lstrm<=length(type))))=[];
%                 d(:,lstrm)=[];
%                 p.link(lstrm,:)=[];
%                 p.link.type = type;
%                 
%                 p.link.type = type;
%                 
%                 if(~ismember('source',p.link.Properties.VariableNames) & ...
%                         ismember('ROI',p.link.Properties.VariableNames))
%                     [p.link,idx] = nirs.util.sortrows(p.link,{'ROI','type'});
%                 else
%                     [p.link,idx] = nirs.util.sortrows(p.link,{'source','detector','type'});
%                 end
%                 
%                 data(i).data  = d(:,idx);
%                 data(i).probe = p;


               clear type;
               
               %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
               % Revisions by bdz 25 Sept 2018
               
               % Allocate a new data matrix because there are only two
               % chromophore concentrations, regardless of # of wavelengths
               
               % determine number of wavelengths measured
               link=p.link;

               [uLambda,~,uLambda_idx] = unique(link.type);

               % If its stored as a cell for some reason, fix this
               if(any(iscell(uLambda)))
                    uLambda=str2num(cell2mat(uLambda));
                    link.type=uLambda(uLambda_idx);
               end


               % grab wavelengths
               extSpectra = [uLambda,nirs.media.getspectra( uLambda )]; % [ wv, hbo, hbr, water, fat]
               
               %preallocate output array
               if(~ismember('source',p.link.Properties.VariableNames) & ...
                       ismember('ROI',p.link.Properties.VariableNames))
                    nchan = unique([p.link.ROI]);
               else
                   nchan = unique([p.link.source p.link.detector]);
               end
               
               % d_chr will store delta chromophore data (d conc)
               % initial size of d_chr will be t x 2 (hbo,hbr) x chan
               d_chr = nan(size(d,1), 2, length(nchan));
               
               % initial size of type_chr will likewise be 2 x chan
               type_chr = cell(2,length(nchan));
               
               % Preserve original index
               link.data_idx=[1:height(link)]';

               % store optode numbers
               link.opt_idx=idx;

               % store type index
               link.uLambda_idx=uLambda_idx;

               % bind spectra
               link.E = extSpectra(link.uLambda_idx,[2,3]); % [ wv, hbo, hbr, water, fat]

               % bind distances
               link.dist=[p.distances];

               % avoid issues with the short (0) seperation values
               %    set to 1mm
               link.dist=max(link.dist,ones(size(link.dist)));

               % drop non-active or other channels
               link=link(~isnan(link.type)&link.type>0,:);
            

               numOptodes=max(idx);
             
                for j = 1:numOptodes
                    % For each optode j

                    link_opt=link(link.opt_idx==j,:);
                    assert( height(link_opt) > 1 )

                    % get data index
                    d_idx=link_opt.data_idx;
                    
                    % wavelengths
                    lambda=link_opt.type;
                  
                    % extinction coefficients
                    E = link_opt.E;
                    
                    % distances
                    L = link_opt.dist;
                    
                    if(isa(obj.PPF,'function_handle'))
                        PPF = obj.PPF(lambda,data(i));      
                    elseif(length(obj.PPF)==1)
                        PPF=repmat(obj.PPF,length(E),1);
                    else
                        PPF=obj.PPF(:);
                    end
                    
                    if(length(E)>2)
                        r=mad(d(:,d_idx)',1,2);
                        r=r-mean(r);
                        s = mad(r, 0) / 0.6745;
                        r = r/s/4.685;
                        w = diag((abs(r)<1) .* (1 - r.^2).^2);
                    else
                        w=diag([1 1]);
                    end
                    
                    % mbll model
                    EL = bsxfun( @times, E, w*L.*PPF );
                    
                    iEL = pinv(EL);
                    
                    % calculates chromophore concentration (uM)
                    d_chr(:,:,j) = (d(:,d_idx)*w*iEL') * 1e6;
                    
                    % new channel type
                    type_chr(:,j) = {'hbo','hbr'};
                    
                end
                
                % p.link needs to be redefined because there are only 2
                % observations per optode-pair, not n-many (for n Lambdas)
                link=p.link;
                types=unique(link.type); 
                link=link(ismember(link.type,types(1:size(type_chr,1))),:);

                
                 p.link = link;

                 l=reshape(type_chr,size(type_chr,1)*size(type_chr,2),1);
                 for idi=1:length(l); ii(idi)=isempty(l{idi}); end;
                l={l{find(~ii)}}';

                p.link.type = l;
                
                if(~ismember('source',p.link.Properties.VariableNames) & ...
                        ismember('ROI',p.link.Properties.VariableNames))
                    [p.link,idx] = nirs.util.sortrows(p.link,{'ROI','type'});
                    p.link=p.link(:,{'ROI','type'});
                else
                    [p.link,idx] = nirs.util.sortrows(p.link,{'source','detector','type'});
                    if(ismember('ShortSeperation',link.Properties.VariableNames))
                        p.link=p.link(idx,{'source','detector','type','ShortSeperation'});
                    else
                        p.link=p.link(idx,{'source','detector','type'});
                    end
                end
                
                % concatenate the channels (3rd dim) into the columns (2nd dim)
                % data is now: t x (species x channel)
                data(i).data  = reshape(d_chr,size(d_chr,1),size(d_chr,2)*size(d_chr,3));
                data(i).probe = p;
                
                % End revisions
                
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
                        extSpectra = nirs.media.getspectra( lambda );
                        
                        clist = [1 2]; % hbo and hbr; need to fix this
                        
                        % extinction coefficients
                        E = extSpectra(:,clist);
                        
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

