classdef bootstrap_result
    
    properties
        truth;
    end
    properties(Dependent=true)
        pvalue;
    end

    properties(Hidden=true)
        ecdf;
        value_bins;
    end

    methods 
        
        function pvalue = get.pvalue(obj)
            pvalue = nan(size(obj.truth));
            for i=1:length(obj.truth)
                if(~isnan(obj.truth(i)))
                    [~,j]=min(abs(obj.value_bins(:,i)-obj.truth(i)));
                    pvalue(i)=2*obj.ecdf(j,i);
                end
            end
            

        end
    end

end
