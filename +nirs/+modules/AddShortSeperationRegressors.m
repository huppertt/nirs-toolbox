classdef AddShortSeperationRegressors < nirs.modules.AbstractModule
    %% AddShortSeperationRegressors - Adds short seperation data as regressors to the GLM model
    %
    
    properties
        baselinecorrect;
        scICA;  % use single channel ICA instead of PCA for defining regressors
    end
    
    methods
        function obj = AddShortSeperationRegressors( prevJob )
            obj.name = 'AddShortSeperationRegressors';
            obj.baselinecorrect=false;
            obj.scICA = false;
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function data = runThis( obj, data )
            for i = 1:numel(data)
                if(~nirs.util.hasshortdistances(data(i)))
                    continue;
                end
                lstss=find( data(i).probe.link.ShortSeperation);
                dd=data(i).data(:,lstss);
                
                if(obj.baselinecorrect)
                    tmp=nirs.core.Data;
                    tmp.time=data(i).time;
                    tmp.data=dd;
                    job=nirs.modules.BaselineCorrection;
                    tmp=job.run(tmp);
                    dd=tmp.data;
                end
                
               dd=dd-ones(size(dd,1),1)*mean(dd,1);
               if(~obj.scICA)
                    dd=orth(dd);
               else
                   dd2=[]; 
                   n=round(10*data(i).Fs);
                   for id=1:size(dd,2); 
                       a=convmtx(dd(:,id),n); 
                       dd2=[dd2 a(:,1:10:end)]; 
                   end;
                   [u,s,v]=svd(dd2);
                   ncom=min(find(cumsum(diag(s))/sum(diag(s))>.8));
                   dd=fastica(dd2', 'approach', 'symm', 'g', 'tanh','numOfIC', ncom)';
                   dd=orth(dd);
               end
                        
                for j=1:size(dd,2)
                    st=nirs.design.StimulusVector;
                    st.name=['SS_PCA' num2str(j)];
                    st.time=data(i).time;
                    st.vector=dd(:,j);
                    st.vector=st.vector-mean(st.vector);
                    st.vector=st.vector./sqrt(var(st.vector));
                    data(i).stimulus(st.name)=st;  
                end
                
            end
        end
    end
    
end

