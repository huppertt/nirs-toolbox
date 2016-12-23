classdef Hyperscanning < nirs.modules.AbstractModule
    %% Hyper-scanning - Computes all-to-all connectivity model between two seperate files
    % Outputs nirs.core.ConnStats object
    
    properties
        corrfcn;  % function to use to compute correlation (see +nirs/+sFC for options)
        divide_events;  % if true will parse into multiple conditions
        min_event_duration;  % minimum duration of events
        link;
        symetric;
    end
    methods
        function obj = Hyperscanning( prevJob )
            obj.name = 'Hypercanning';
            obj.corrfcn = @(data)nirs.sFC.ar_corr(data,'4xFs',true);  %default to use AR-whitened robust correlation
            obj.divide_events=false;
            obj.min_event_duration=30;
            obj.symetric=true;
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function connStats = runThis( obj, data )
            
            
            if(isempty(obj.link))
                % NIRx files have a hyperscan variable upon loading that I
                % can use here
                if(ismember('hyperscan',nirs.createDemographicsTable(data).Properties.VariableNames))
                    hyperscanfiles=nirs.createDemographicsTable(data).hyperscan;
                    
                    for i=1:length(hyperscanfiles); if(isempty(hyperscanfiles{i})); hyperscanfiles{i}=''; end; end;
                    uniquefiles=unique(cellstr(vertcat(hyperscanfiles{:})));
                    [ia,ib]=ismember(hyperscanfiles,uniquefiles);
                    
                    for i=1:length(uniquefiles)
                        lst=find(ib==i);
                        ScanA(i,1)=lst(1);
                        ScanB(i,1)=lst(2);
                    end
                    
                    OffsetA = zeros(size(ScanA));  % The time shift of the "A" files (in sec)
                    OffsetB = zeros(size(ScanB));  % The time shift of the "B" files (in sec)
                    
                    link = table(ScanA,ScanB,OffsetA,OffsetB);
                    obj.link=link;
                else
                    warning('link variable must be specified');
                    connStats=nirs.core.sFCStats();
                    return
                end
            end
            
            
            
            for i=1:height(obj.link)
                idxA = obj.link.ScanA(i);
                idxB = obj.link.ScanB(i);
                
                dataA = data(idxA).data;
                timeA = data(idxA).time+obj.link.OffsetA(i);
                dataB = data(idxB).data;
                timeB = data(idxB).time+obj.link.OffsetB(i);
                
                % Make sure we are using the same time base
                time=[max(timeA(1),timeB(1)):1/data(idxA).Fs:min(timeA(end),timeB(end))];
                for id=1:size(dataA,2)
                    dataA(1:length(time),id)=interp1(timeA,dataA(:,id),time);
                    dataB(1:length(time),id)=interp1(timeB,dataB(:,id),time);
                end
                dataA=dataA(1:length(time),:);
                dataB=dataB(1:length(time),:);
                
                dataA=dataA-ones(length(time),1)*mean(dataA,1);
                dataB=dataB-ones(length(time),1)*mean(dataB,1);
                
                dataA=dataA./(ones(length(time),1)*mad(dataA,1,1));
                dataB=dataB./(ones(length(time),1)*mad(dataB,1,1));
                
                connStats(i)=nirs.core.sFCStats();
                connStats(i).type = obj.corrfcn;
                connStats(i).description= data(i).description;
                connStats(i).probe=data(idxA).probe;
                connStats(i).demographics={data(idxA).demographics; data(idxB).demographics};
                
                cond={};
                if(obj.divide_events)
                    stimNames=unique(horzcat(data(idxA).stimulus.keys, data(idxB).stimulus.keys));
                    stim=Dictionary;
                    for idx=1:length(stimNames)
                        onsets=[];
                        dur=[];
                        s=data(idxA).stimulus(stimNames{idx});
                        if(~isempty(s))
                            onsets=[onsets s.onset];
                            dur=[dur s.dur];
                        end
%                         s=data(idxB).stimulus(stimNames{idx});
%                         if(~isempty(s))
%                             onsets=[onsets; s.onset];
%                             dur=[dur; s.dur];
%                         end
                        ss=nirs.design.StimulusEvents;
                        ss.name=stimNames{idx};
                        ss.onset=onsets;
                        ss.dur=dur;
                        ss.amp=ones(size(dur));
                        stim(stimNames{idx})=ss;
                    end
                    
                    
                    Basis=Dictionary;
                    Basis('default')=nirs.design.basis.BoxCar;
                    [X, names] = nirs.design.createDesignMatrix(stim,time,Basis);
                    mask={};
                    if(~isempty(find(sum(abs(X),2)==0)))
                        mask{1}=1*(sum(abs(X),2)==0);
                        cond{1}='rest';
                    end
                    
                    for idx=1:size(X,2);
                        if(~all(X(:,idx)==0))
                            mask{end+1}=(X(:,idx)>0)*1;
                            cond{end+1}=names{idx};
                        else
                            disp(['discluding stimulus: ' names{idx}]);
                        end
                    end
                    
                else
                    mask={ones(length(time),1)};
                    cond={'rest'};
                end
                
                connStats(i).R=[];
                for cIdx=1:length(mask)
                    tmp=data(idxA);
                    
                    dd=[dataA dataB]+(1i)*(mask{cIdx}*ones(1,size(dataA,2)*2));
                    if(obj.symetric)
                        Z=zeros(40,size(dd,2));
                        dd=[dd; Z; [dataB dataA]+(1i)*(mask{cIdx}*ones(1,size(dataA,2)*2))];
                    end
                    
                    [ii,~]=find(dd~=0);
                    lstBad=[max(ii)+1:length(time)];
                    lstBad=[lstBad 1:min(ii)]; 
                    dd(lstBad,:)=[];
                    
                    tmp.time=time;
                    tmp.time(lstBad)=[];
                    tmp.data=dd;
                   
                    [r,p,dfe]=obj.corrfcn(tmp);
                    
                    if(cIdx==1)
                        connStats(i).probe = nirs.util.createhyperscanprobe(connStats(i).probe);
                    end
                    r(1:end/2,1:end/2)=0;
                    r(1+end/2:end,1+end/2:end)=0;
                    
                    
                    connStats(i).dfe(cIdx)=dfe;
                    connStats(i).R(:,:,cIdx)=r;
                    
                end
                connStats(i).conditions=cond;
                disp(['Finished ' num2str(i) ' of ' num2str(height(obj.link))])
                
                
            end
            
        end
    end
end
        
