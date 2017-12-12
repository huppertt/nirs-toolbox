classdef AddShortSeperationRegressors < nirs.modules.AbstractModule
    %% AddShortSeperationRegressors - Adds short seperation data as regressors to the GLM model
    %
    
    properties
          scICA;  % use single channel ICA instead of PCA for defining regressors
    end
    
    methods
        function obj = AddShortSeperationRegressors( prevJob )
            obj.name = 'AddShortSeperationRegressors';
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
                
                dd=dd-ones(size(dd,1),1)*mean(dd,1);
               
                if(~obj.scICA)
                    dd=orth(dd);
                else
                   dd2=[dd, [diff(dd); zeros(1,size(dd,2))],[diff(diff(dd)); zeros(2,size(dd,2))]]; 
                   
                   dd=orth(dd2);
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

