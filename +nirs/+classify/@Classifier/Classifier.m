classdef Classifier
    %% CLASSIFIER - holds data for NIRS/EEG/multimodal classifier models
    %
    
    properties
        model;
        training_data;
        training_labels
    end
    
    properties(Dependent = true)
        accuracy;
    end
    properties( Hidden=true)
        trainedmodel;
    end
    
    
    methods 
        function obj = Classifier(varargin)
            
            if(nargin>0 && (isa(varargin{1},'cell') | isa(varargin{1},'eeg.core.ChannelStats') | isa(varargin{1},'nirs.core.ChannelStats')))
                [training_data,training_labels]=nirs.classify.util.stats2labels(varargin{1});
                if(nargin>1)
                    model=varargin{2};
                else
                    model='treebagger';
                end
            
            else
                
                if(nargin<1)
                    training_data=[];
                else
                    training_data=varargin{1};
                end
                if(nargin<2)
                    training_labels={};
                else
                    training_labels=varargin{2}
                end
                if(nargin<3 || isempty(varargin{3}))
                    model='treebagger';
                else
                    model=varargin{3};
                end
            end
            
            
            obj.model=model;
            obj.training_data=training_data;
            obj.training_labels=training_labels;
            obj=trainmodel(obj);
        end
        
        function obj = set.model( obj, model )
            allowed={'treebagger','svm','lda'};
            
            oldmodel=obj.model;
            
          if(~ ismember(lower(model),allowed) )
              disp('Model must be one of the following:')
              disp(allowed)
              return
          end
           obj.model = model;
           if(~strcmp(lower(oldmodel),lower(model)))
                obj=trainmodel(obj);
           end
        end
        
        function obj = trainmodel(obj)
            if(isempty(obj.training_data) || isempty(obj.training_labels))
                return
            end
            switch(lower(obj.model))
                case('treebagger')
                    n=length(unique(obj.training_labels))*2;
                    obj.trainedmodel=TreeBagger(n,obj.training_data,obj.training_labels,'OOBPrediction','on');
            end
            
        end
        
        function acc = get.accuracy(obj)
            pred=obj.predict(obj.training_data);
            if(isa(pred,'cell'))
                acc=length(find(strcmp(pred,obj.training_labels)))/length(pred);
            else
                acc=length(find(pred==obj.training_labels))/length(pred);
            end        
        end
        
        function pred = predict(obj,data)
            
            if( (isa(data,'cell') && (isa(data{1},'nirs.core.ChannelStats') || isa(data{1},'eeg.core.ChannelStats'))) ||...
                     isa(data,'nirs.core.ChannelStats') || isa(data,'eeg.core.ChannelStats') )
                data=nirs.classify.util.stats2labels(data);
            elseif( (isa(data,'cell') && (isa(data{1},'nirs.core.Data') || isa(data{1},'eeg.core.Data'))) ||...
                    isa(data,'nirs.core.Data') || isa(data,'eeg.core.Data') )
                if(isa(data,'cell'))
                    for i=1:length(data)
                        d{i}=vertcat(data{i}(:).data);
                    end
                    data=horzcat(d{:});
                else
                    data=vertcat(data(:).data);
                end
            end
                
            
            if(~isempty(obj.trainedmodel))
                pred=obj.trainedmodel.predict(data);
            else
                pred=[];
            end
                    
        end
        
    end
end