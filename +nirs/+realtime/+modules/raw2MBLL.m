classdef raw2MBLL < nirs.realtime.modules.AbstractModule
    % real-time implementation of optical density to HbX
    properties
       PPF=6;
    end
    
    properties(Hidden=true)
       mbll=[];
    end
    
    methods
        function obj=raw2MBLL(prevJob)
            obj.name='RT-MBLL';
            if nargin > 0
                obj.prevJob = prevJob;
            end
            obj.hpf.lowpass=[];
            obj.hpf.highpass=0.001;
        end
        
        function obj=resetThis(obj)
            obj.mbll=[];
        end
        
        function [hb,t,probe,stimulus] = updateThis(obj,d,t,probe,stimulus)
            
            if(isempty(obj.mbll))
               lambda=unique(probe.link.type);
               [~,~,idx] = nirs.util.uniquerows(table([probe.link.source probe.link.detector]));
                
                link=table;
               
                for j = 1:max(idx)
                    lst = idx == j;
                    lst=find(lst);
                    link=[link; probe.link(lst(1:2),:)];
                    assert( length(lst) > 1 )
                    lambda = probe.link.type(lst);
                   
                    ext = nirs.media.getspectra( lambda );
                    clist = [1 2]; % hbo and hbr; need to fix this
                    
                    % extinction coefficients
                    E = ext(:,clist);
                    
                    % distances
                    L = probe.distances(lst);
                    L=max(L,1);  % avoid issues with the short (0) seperation values
                    
                    if(isa(obj.PPF,'function_handle'))
                        PPF = obj.PPF(lambda,data(i));      
                    elseif(length(obj.PPF)==1)
                        PPF=repmat(obj.PPF,length(lambda),1);
                    else
                        PPF=obj.PPF(:);
                    end
                    
                    if(length(lambda)>2)
                        r=mad(d(:,lst)',1,2);
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
                    obj.mbll.A{j}=w*iEL' * 1e6;
                    obj.mbll.lst{j}=lst;
                    obj.mbll.lst2{j}=[2*(j-1)+1 2*(j-1)+2];
                    % new channel type
                    type_chr(:,j) = {'hbo','hbr'};
                    
                    obj.mbll.cnt=1;
                    obj.mbll.dc=d;
                    
                end
                probe.link = link;
                probe.link.type = reshape(type_chr,size(type_chr,1)*size(type_chr,2),1);
                [probe.link,idx] = nirs.util.sortrows(probe.link,{'source','detector','type'});
                obj.mbll.probe=probe;
                
            end
            probe=obj.mbll.probe;
            hb=zeros(size(d,1),height(obj.mbll.probe.link));
            
            
            for jj=1:size(d,1)
                obj.mbll.dc=(obj.mbll.dc*obj.mbll.cnt+d(jj,:))/(obj.mbll.cnt+1);
                obj.mbll.cnt=obj.mbll.cnt+1;
                d(jj,:)=d(jj,:)./obj.mbll.dc;
            end
            
            d=-log(d);
            

            
            
            for id=1:length(obj.mbll.lst)
                hb(:,obj.mbll.lst2{id})=d(:,obj.mbll.lst{id})*obj.mbll.A{id};
            end
            
            
        end
            
    end
    
end

