classdef OrdinalVector
    %This class holds ordinal vectors
    
    properties
        name
        numericvalues;
        classlabels;
        classedges;
        time
    end
    
    methods
        function vec = getStimVector( obj, time )
            vec = interp1( obj.time, obj.numericvalues, time );
            
            if(~isempty(obj.classedges))
                vec = ordinal(vec,obj.classlabels,[],obj.classedges);
            else
                vec = ordinal(vec,obj.classlabels);
            end
            
        end
        
        function draw(obj)
            if(~isempty(obj.classedges))
                vec = ordinal(obj.numericvalues,obj.classlabels,[],obj.classedges);
            else
                vec = ordinal(obj.numericvalues,obj.classlabels);
            end
            figure;
            plot(obj.time,vec);
            xlabel('Time (sec)')
            ylabel(obj.name);
            set(gca,'YTick',[1:length(obj.classlabels)]);
            set(gca,'YTickLabel',obj.classlabels);
            
        end
        
        function st=convert2categorical(obj)
            st = nirs.design.CategoricalVector;
            st.name=obj.name;
            st.values = categorical(ordinal(obj.numericvalues,obj.classlabels,[],obj.classedges));
            st.time=obj.time;
        end
            
        
        
    end
    
end