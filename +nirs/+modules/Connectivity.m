classdef Connectivity < nirs.modules.AbstractModule
%% CONNECTIVITY - Computes all-to-all connectivity model.
% Outputs nirs.core.ConnStats object

    properties
        corrfcn;  % function to use to compute correlation (see +nirs/+sFC for options)
        divide_events;  % if true will parse into multiple conditions
        min_event_duration;  % minimum duration of events
        ignore;
    end
    methods
        function obj = Connectivity( prevJob )
           obj.name = 'Connectivity';
           obj.corrfcn = @(data)nirs.sFC.ar_corr(data,'4xFs',true);  %default to use AR-whitened robust correlation
           obj.divide_events=false;
           obj.min_event_duration=30;
           obj.ignore=10;
           if nargin > 0
               obj.prevJob = prevJob;
           end
           obj.citation='Santosa, H., Aarabi, A., Perlman, S. B., & Huppert, T. J. (2017). Characterization and correction of the false-discovery rates in resting state connectivity using functional near-infrared spectroscopy. Journal of Biomedical Optics, 22(5), 055002-055002.';
           
        end
        
        function connStats = runThis( obj, data )
            for i = 1:numel(data)
                
               
                connStats(i)=nirs.core.sFCStats();
                connStats(i).type = obj.corrfcn;
                connStats(i).description= ['Connectivity model of ' data(i).description];
                connStats(i).probe=data(i).probe;
                connStats(i).demographics=data(i).demographics;
     
                mask={}; cond={};
                if(obj.divide_events)
                    
                    stim=data(i).stimulus;
                    cnt=1;
                    for idx=1:length(stim.keys)
                        s=stim(stim.keys{idx});
                        lst=find(s.dur-2*obj.ignore>obj.min_event_duration);
                        if(length(lst)>0)
                            s.onset=s.onset(lst);
                            s.dur=s.dur(lst);
                            
                            disp(['Spliting condition: ' stim.keys{idx}]);
                            r=zeros(size(data(i).data,2),size(data(i).data,2),length(s.onset));
                            dfe=zeros(length(s.onset),1);
                            for j=1:length(s.onset)
                                disp(['   ' num2str(j) ' of ' num2str(length(s.onset))]);
                                lstpts=find(data(i).time>s.onset(j)+obj.ignore &...
                                    data(i).time<s.onset(j)+s.dur(j)-obj.ignore);
                                tmp=data(i);
                                tmp.data=data(i).data(lstpts,:);
                                tmp.time=data(i).time(lstpts);
                                [r(:,:,j),p,dfe(j)]=obj.corrfcn(tmp);
                            end
                            
                            connStats(i).dfe(cnt)=sum(dfe);
                            connStats(i).R(:,:,cnt)=tanh(mean(atanh(r),3));
                            connStats(i).conditions{cnt}=stim.keys{idx};
                            cnt=cnt+1;
                        else
                            disp(['Skipping condition: ' stim.keys{idx} ...
                                ': No events > ' num2str(2*obj.ignore+obj.min_event_duration) 's']);
                        end
                        
                    end

                else
                    tmp=data(i);
                    lst=find(tmp.time<obj.ignore | tmp.time>tmp.time(end)-obj.ignore);
                    tmp.data(lst,:)=[];
                    tmp.time(lst)=[];
                    [r,p,dfe]=obj.corrfcn(tmp);
                    
                    connStats(i).dfe=dfe;
                    connStats(i).R=r;
                    connStats(i).conditions=cellstr('rest');
                end
                
                
                disp(['Finished ' num2str(i) ' of ' num2str(length(data))]);
            end
        end
        
    end
    
end
