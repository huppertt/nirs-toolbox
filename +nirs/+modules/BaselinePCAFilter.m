classdef BaselinePCAFilter < nirs.modules.AbstractModule
    %% PCA filter from rest data
    
    properties
        nSV;  % Number of SV to remove
        table;  %
        splittypes;
        discard;
    end
    
    methods
        
        function obj = BaselinePCAFilter( prevJob )
            obj.name = 'Baseline PCA filter';
            obj.nSV = 1;
            obj.discard =  true;
            splittypes= true;
            obj.table =table(1,1,'VariableNames',{'Task','Baseline'});
            obj.table(1,:)=[];
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function datafilt = runThis( obj, data )
            uniqBaseline = unique(obj.table.Baseline);
            datafilt=data;
            
            types=unique(data(1).probe.link.type);
            for tI=1:length(types)
                
                if(tI>1 & ~obj.splittypes)
                    continue;
                end
                
                if(obj.splittypes)
                    lst=find(ismember(data(1).probe.link.type,types{tI}));
                else
                    lst=1:size(data(1).data,2);
                end
                
                for idx=1:length(uniqBaseline)
                    thesescans = obj.table.Task(find(obj.table.Baseline==uniqBaseline(idx)));
                    
                    % [SVD]= PCA(data(uniqBaseline(idx))) << compute PCA from
                    
                    dBase=data(uniqBaseline(idx)).data(:,lst);
                    for ii=1:size(dBase,2)
                        dBase(:,ii)=detrend(dBase(:,ii));
                    end
                    
                    % svd
                    [u, s, v] = svd(dBase'*dBase,'econ');
                    %s = diag(s);
                    
                    if(obj.nSV<1)
                        norms=cumsum(diag(s))/sum(diag(s));
                        nSV=max(find(norms<obj.nSV));
                        lstncomp=1:nSV;
                        if(~isempty(nSV))
                            disp(['Removing ' num2str(nSV) ' of ' num2str(length(diag(s))) ' (' num2str(norms(nSV)*100) '%)']);
                        else
                            disp('not removing any components');
                        end
                    else
                        lstncomp = 1:obj.nSV;
                    end
                    % baseline
                    for idx2=1:length(thesescans)
                        % datafiltered(thesescans(idx2))=PCAfilt(SVD,data(thesescans(idx2)))
                        dTask= data(thesescans(idx2)).data(:,lst);
                        dTask=detrend(dTask);
                        u = dTask*v*inv(s);
                        
                        datafilt(thesescans(idx2)).data(:,lst) = dTask - u(:,lstncomp)*s(lstncomp,lstncomp)*v(:,lstncomp)';
                        
                    end
                end
            end
            if(obj.discard)
                % remove non-filtered scans
                lst=[1:length(datafilt)];
                datafilt(~ismember(lst,obj.table.Task))=[];
            end
            
        end
    end
    
end

