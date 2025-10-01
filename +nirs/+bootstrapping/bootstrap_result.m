classdef bootstrap_result
    
    properties
        truth;
    end
    properties(Dependent=true)
        pvalue;
        pvalue_lower_bound;
        pvalue_upper_bound;
    end

    properties(Hidden=true)
        pcdf_estimate;
        pcdf_lower_bounds;
        pcdf_upper_bounds;
        value_bins;
    end

    methods 
        function pvalue_lower_bound = get.pvalue_lower_bound(obj)
            LB = zeros(size(obj.truth));
            UB = zeros(size(obj.truth));
            
            for i=1:size(obj.truth,1)
                for j=1:size(obj.truth,2)
                    idx1=min(find(obj.value_bins(:,i,j)>=obj.truth(i,j)));
                    idx2=max(find(obj.value_bins(:,i,j)<=obj.truth(i,j)));
                    lst=[idx1 idx2];
                    
                    LB(i,j)=nanmin(obj.pcdf_lower_bounds(lst,i,j));
                    UB(i,j)=nanmax(obj.pcdf_upper_bounds(lst,i,j));

                    if(obj.truth(i,j)>0)
                        UB(i,j)=1-UB(i,j);
                        LB(i,j)=1-LB(i,j);
                    end
                end
            end
                
            pvalue_lower_bound=nanmin(UB,LB);
            pvalue_upper_bound=nanmax(UB,LB);
        end
        function pvalue_upper_bound = get.pvalue_upper_bound(obj)
            LB = zeros(size(obj.truth));
            UB = zeros(size(obj.truth));
            
            for i=1:size(obj.truth,1)
                for j=1:size(obj.truth,2)
                    idx1=min(find(obj.value_bins(:,i,j)>=obj.truth(i,j)));
                    idx2=max(find(obj.value_bins(:,i,j)<=obj.truth(i,j)));
                    lst=[idx1 idx2];
                    
                    LB(i,j)=nanmin(obj.pcdf_lower_bounds(lst,i,j));
                    UB(i,j)=nanmax(obj.pcdf_upper_bounds(lst,i,j));

                    if(obj.truth(i,j)>0)
                        UB(i,j)=1-UB(i,j);
                        LB(i,j)=1-LB(i,j);
                    end
                end
            end
                
            pvalue_lower_bound=nanmin(UB,LB);
            pvalue_upper_bound=nanmax(UB,LB);
        end
        function pvalue = get.pvalue(obj)
            pvalue = zeros(size(obj.truth));
            for i=1:size(obj.truth,1)
                for j=1:size(obj.truth,2)
                    idx1=min(find(obj.value_bins(:,i,j)>=obj.truth(i,j)));
                    idx2=max(find(obj.value_bins(:,i,j)<=obj.truth(i,j)));
                    lst=[idx1 idx2];
                    pvalue(i,j)=nanmean(obj.pcdf_estimate(lst,i,j));
                     if(obj.truth(i,j)<0)
                         pvalue(i,j)=1-pvalue(i,j);
                     end

                end
            end
            pvalue=1-pvalue;
        end
    end

end
