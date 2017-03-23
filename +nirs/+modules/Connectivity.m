classdef Connectivity < nirs.modules.AbstractModule
%% CONNECTIVITY - Computes all-to-all connectivity model.
% Outputs nirs.core.ConnStats object

    properties
        corrfcn;  % function to use to compute correlation (see +nirs/+sFC for options)
        divide_events;  % if true will parse into multiple conditions
        min_event_duration;  % minimum duration of events
    end
    methods
        function obj = Connectivity( prevJob )
           obj.name = 'Connectivity';
           obj.corrfcn = @(data)nirs.sFC.ar_corr(data,'4xFs',true);  %default to use AR-whitened robust correlation
           obj.divide_events=false;
           obj.min_event_duration=30;
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function connStats = runThis( obj, data )
            for i = 1:length(data)
                
               
                connStats(i)=nirs.core.sFCStats();
                connStats(i).type = obj.corrfcn;
                connStats(i).description= ['Connectivity model of ' data(i).description];
                connStats(i).probe=data(i).probe;
                connStats(i).demographics=data(i).demographics;
     
                mask={}; cond={};
                if(obj.divide_events)
                    
                    stim=data(i).stimulus;
                    for idx=1:length(stim.keys)
                        s=stim(stim.keys{idx});
                        lst=find(s.dur<obj.min_event_duration);
                        s.onset(lst)=[];
                        s.dur(lst)=[];
                        s.amp(lst)=[];
                        s.amp(:)=1;
                        stim(stim.keys{idx})=s;
                    end
                    
                    Basis=Dictionary;
                    Basis('default')=nirs.design.basis.BoxCar;
                    [X, names] = nirs.design.createDesignMatrix(stim,data(i).time,Basis);
                    mask{1}=1*(sum(abs(X),2)==0);
                    cond{1}='rest';
                    
                    for idx=1:size(X,2);
                        if(~all(X(:,idx)==0))
                            mask{end+1}=X(:,idx);
                            cond{end+1}=names{idx};
                        else
                            disp(['discluding stimulus: ' names{idx}]);
                        end
                    end
                
                else
                    mask={ones(size(data(i).data,1),1)};
                    cond={'rest'};
                end
                
                for cIdx=1:length(mask)
                    tmp=data(i);
                    tmp.data=tmp.data-ones(size(tmp.data,1),1)*mean(tmp.data,1);
                    tmp.data=tmp.data-ones(size(tmp.data,1),1)*mean(tmp.data,1);
                    tmp.data=tmp.data+(1i)*(mask{cIdx}*ones(1,size(tmp.data,2)));
                    [r,p,dfe]=obj.corrfcn(tmp);
                    
                    connStats(i).dfe(cIdx)=dfe;
                    connStats(i).R(:,:,cIdx)=r;
                    
                end
                connStats(i).conditions=cond;
                disp(['Finished ' num2str(i) ' of ' num2str(length(data))]);
            end
        end
        
    end
    
end
