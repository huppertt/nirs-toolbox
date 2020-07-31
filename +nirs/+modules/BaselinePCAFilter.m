classdef BaselinePCAFilter < nirs.modules.AbstractModule
    %% PCA filter from rest data
    
    properties
        nSV;  % Number of SV to remove
        table;  %
        splittypes;
        discard;
        keep_components=false;
    end
    
    methods
        
        function obj = BaselinePCAFilter( prevJob )
            obj.name = 'Baseline PCA filter';
            obj.nSV = 1;
            obj.discard =  true;
            obj.splittypes= true;
            obj.table =table(1,1,'VariableNames',{'Task','Baseline'});
            obj.table(1,:)=[];
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function datafilt = runThis( obj, data )
            uniqBaseline = unique(obj.table.Baseline);
            datafilt=data;
            
            if(obj.keep_components)
                for i=1:length(datafilt)
                    comp=datafilt(i);
                    comp.auxillary('eigenvalues')=[];
                    datafilt(i).auxillary('PCAfilter')=comp;
                end
            end
            
            types=unique(data(1).probe.link.type);
            for tI=1:length(types)
                
                if(tI>1 & ~obj.splittypes)
                    continue;
                end
                
                if(~iscell(types))
                    types=num2cell(types);
                end
                
                if(obj.splittypes)
                    lst=find(ismember(data(1).probe.link.type,types{tI}));
                else
                    lst=1:size(data(1).data,2);
                end
                
                if(~isempty(uniqBaseline))
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
                            
                            m = mean(data(thesescans(idx2)).data(:,lst),1);
                            data(thesescans(idx2)).data(:,lst) = bsxfun(@minus, data(thesescans(idx2)).data(:,lst), m);
                            
                            
                            % datafiltered(thesescans(idx2))=PCAfilt(SVD,data(thesescans(idx2)))
                            dTask= data(thesescans(idx2)).data(:,lst);
                            dTask=detrend(dTask);
                            u = dTask*v*inv(s);
                            
                            datafilt(thesescans(idx2)).data(:,lst) = dTask - u(:,lstncomp)*s(lstncomp,lstncomp)*v(:,lstncomp)';
                            datafilt(thesescans(idx2)).data(:,lst) = bsxfun(@plus, datafilt(thesescans(idx2)).data(:,lst), m);
                            
                            if(obj.keep_components)
                                comp=datafilt(thesescans(idx2)).auxillary('PCAfilter');
                                comp.data(:,lst) = u(:,lstncomp)*s(lstncomp,lstncomp)*v(:,lstncomp)';
                                eig=comp.auxillary('eigenvalues');
                                eig(:,end+1)=diag(s);
                                comp.auxillary('eigenvalues')=eig;
                                datafilt(thesescans(idx2)).auxillary('PCAfilter')=comp;
                            end
                            
                        end
                    end
                    
                    
                else
                    %use the baseline rest period from the data file itself
                    for i=1:numel(data)
                        m = mean(data(i).data(:,lst),1);
                        data(i).data(:,lst) = bsxfun(@minus, data(i).data(:,lst), m);
                        
                        basis=nirs.design.basis.BoxCar;
                        basis.lagTime=0;
                        basis.irf_dur=12;
                        b=Dictionary;
                        b('default')=basis;
                        X=nirs.design.createDesignMatrix(data(i).stimulus,data(i).time,b);
                        lsts=find(sum(abs(X),2)==0);
                        dlst=[1; find(diff([lsts; size(X,1)+1])>1); length(lsts)+1];
                        dBase=data(i).data(lsts,lst);
                        for j=1:length(dlst)-1
                            dBase(dlst(j):dlst(j+1)-1,:)= dBase(dlst(j):dlst(j+1)-1,:)-ones(dlst(j+1)-dlst(j),1)*mean(dBase(dlst(j):dlst(j+1)-1,:),1);
                        end
                        
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
                        dTask= data(i).data(:,lst);
                        
                        dTask=detrend(dTask);
                        u = dTask*v*pinv(s);
                        
                        datafilt(i).data(:,lst) = dTask - u(:,lstncomp)*s(lstncomp,lstncomp)*v(:,lstncomp)';
                        datafilt(i).data(:,lst)= bsxfun(@plus, datafilt(i).data(:,lst), m);
                        
                        
                        if(obj.keep_components)
                            comp=datafilt(i).auxillary('PCAfilter');
                            comp.data(:,lst) = u(:,lstncomp)*s(lstncomp,lstncomp)*v(:,lstncomp)';
                            eig=comp.auxillary('eigenvalues');
                            eig(:,end+1)=diag(s);
                            comp.auxillary('eigenvalues')=eig;
                            datafilt(i).auxillary('PCAfilter')=comp;
                        end
                        
                    end
                    
                    
                end
            end
            if(obj.discard & (~isempty(uniqBaseline)))
                % remove non-filtered scans
                lst=[1:length(datafilt)];
                datafilt(~ismember(lst,obj.table.Task))=[];
            end
        end
    end
end

