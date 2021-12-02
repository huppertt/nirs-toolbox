classdef AddAuxRegressors < nirs.modules.AbstractModule
    %% AddAuxRegressors - Adds auxillary data as regressors to the GLM model
    %
    
    properties
        normalize;  % normalize the regressors
        orth;
        label;
    end
    
    methods
        function obj = AddAuxRegressors( prevJob )
            obj.name = 'AddAuxRegressors';
            obj.normalize = true;
            obj.label={};
            obj.orth = true;
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function data = runThis( obj, data )
            for i = 1:numel(data)
                
                auxlabels=data(i).auxillary.keys;
                idx=find(ismember(lower(auxlabels),lower(obj.label)));
                if(ismember('all',lower(obj.label)))
                    idx=1:data(i).auxillary.count;
                end
                for j=1:length(obj.label)
                    if(~isempty(strfind(obj.label{j},'*')))
                        ll=strsplit(obj.label{j},'*');
                        for k=1:length(auxlabels)
                            found=true;
                            for l=1:length(ll)
                                if(~isempty(ll{l}))
                                    found=found & ~isempty(strfind(auxlabels{k},ll{l}));
                                end
                            end
                            if(found)
                                idx=[idx k];
                            end
                        end
                    end
                end
                
                if(isempty(idx))
                    continue;
                end
                if(obj.orth)
                    dd=[];
                    for j=1:length(idx)
                        
                        aux=data(i).auxillary(data(i).auxillary.keys{idx(j)});
                        if(length(aux.time)~=length(unique(aux.time)))
                            [aux.time,ia]=unique(aux.time);
                            aux.data=aux.data(ia,:);
                        end
                        for k=1:size(aux.data,2)
                            
                            
                            
                            dd=[dd interp1(aux.time,aux.data(:,k),data(i).time,'linear','extrap')];
                        end
                    end
                    if(obj.normalize)
                        dd=dd-ones(size(dd,1),1)*nanmean(dd,1);
                        dd=dd./(ones(size(dd,1),1)*nanstd(dd,[],1));
                    end
                    dd=orth(dd);
                    for j=1:size(dd,2)
                        st=nirs.design.StimulusVector;
                        st.regressor_no_interest=true;
                        st.name=['AuxPCA' num2str(j)];
                        st.time=data(i).time;
                        st.vector=dd(:,j);
                        data(i).stimulus(st.name)=st;
                    end
                else
                    for j=1:length(idx)
                        dd=[];
                        aux=data(i).auxillary(data(i).auxillary.keys{idx(j)});
                        for k=1:size(aux.data,2)
                            dd=[dd interp1(aux.time,aux.data(:,k),data(i).time)];
                        end
                        
                        if(obj.normalize)
                            dd=dd-ones(size(dd,1),1)*mean(dd,1);
                            dd=dd./(ones(size(dd,1),1)*std(dd,[],1));
                        end
                        dd=orth(dd);
                        for k=1:size(dd,2)
                            st=nirs.design.StimulusVector;
                            st.regressor_no_interest=true;
                            st.name=[data(i).auxillary.keys{idx(j)} num2str(k)];
                            st.time=data(i).time;
                            st.vector=dd(:,k);
                            data(i).stimulus(st.name)=st;
                        end
                    end
                end
                
                
                
                
            end
        end
    end
    
end

