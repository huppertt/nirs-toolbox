classdef Hyperscanning < nirs.modules.AbstractModule
    %% Hyper-scanning - Computes all-to-all connectivity model between two seperate files
    % Outputs nirs.core.ConnStats object
    
    properties
        corrfcn;  % function to use to compute correlation (see +nirs/+sFC for options)
        divide_events;  % if true will parse into multiple conditions
        min_event_duration;  % minimum duration of events
        link;
        symetric;
        verbose;
        ignore;  % seconds at start/end of each scan or block to ignore
        multitrial; % flag to enable computing once per condition rather than averaging over blocks/trials
    end
    properties(Hidden=true)
        % I like to hide options that I don't want the average person
        % playing with
        
        cache_dir;  % (optional) directory to cache results (unset disables caching)
        cache_rebuild;  % (optional) force rebuild of cached results (don't load previous results from cache, only save new results)
    end
    methods
        function obj = Hyperscanning( prevJob )
            obj.name = 'Hyperscanning';
            obj.corrfcn = @(data)nirs.sFC.ar_corr(data,'18xFs',true);  %default to use AR-whitened robust correlation
            obj.divide_events=false;
            obj.min_event_duration=30;
            obj.symetric=true;
            obj.ignore=10;
            obj.multitrial=false; 
            obj.cache_dir='';
            obj.verbose=false;
            obj.cache_rebuild=false;
            if nargin > 0
                obj.prevJob = prevJob;
            end
            
            obj.citation='Santosa, H., Aarabi, A., Perlman, S. B., & Huppert, T. J. (2017). Characterization and correction of the false-discovery rates in resting state connectivity using functional near-infrared spectroscopy. Journal of Biomedical Optics, 22(5), 055002-055002.';
            
        end
        
        function connStats = runThis( obj, data )
            
            
            if(isempty(obj.link))
                % NIRx files have a hyperscan variable upon loading that I
                % can use here
                if(ismember('hyperscan',nirs.createDemographicsTable(data).Properties.VariableNames))
                    hyperscanfiles=nirs.createDemographicsTable(data).hyperscan;
                    
                    for i=1:length(hyperscanfiles); if(isempty(hyperscanfiles{i})); hyperscanfiles{i}=''; end; end;
                    uniquefiles=unique(hyperscanfiles);
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
            
            connStats(1:height(obj.link))=nirs.core.sFCStats();
            
            for i=1:height(obj.link)
                
                idxA = obj.link.ScanA(i);
                idxB = obj.link.ScanB(i);
                
                dataA = data(idxA).data;
                timeA = data(idxA).time+obj.link.OffsetA(i);
                dataB = data(idxB).data;
                timeB = data(idxB).time+obj.link.OffsetB(i);
                
                % Make sure we are using the same time base
                if isequal(timeA,timeB)
                    time = timeA;
                else
                    time=[max(timeA(1),timeB(1)):1/data(idxA).Fs:min(timeA(end),timeB(end))];
                    for id=1:size(dataA,2)
                        dataA(1:length(time),id)=interp1(timeA,dataA(:,id),time);
                        dataB(1:length(time),id)=interp1(timeB,dataB(:,id),time);
                    end
                    dataA=dataA(1:length(time),:);
                    dataB=dataB(1:length(time),:);
                end
                
                connStats(i).type = obj.corrfcn;
                connStats(i).description= data(idxA).description;
                connStats(i).probe = nirs.util.createhyperscanprobe(data(idxA).probe);
                connStats(i).demographics={data(idxA).demographics; data(idxB).demographics};
                connStats(i).R=[];
                
                if(obj.divide_events)
                    stim=data(idxA).stimulus;
                    cnt=1;
                    for idx=1:length(stim.keys)
                        s=stim(stim.keys{idx});
                        lst=find(s.dur-2*obj.ignore>obj.min_event_duration);
                        if(length(lst)>0)
                            s.onset=s.onset(lst);
                            s.dur=s.dur(lst);
                            
                            if(obj.verbose)
                                disp(['Spliting condition: ' stim.keys{idx}]);
                            end
                            n1=size(dataA,2)+size(dataB,2);
                            r=zeros(n1,n1,length(s.onset));
                            dfe=zeros(length(s.onset),1);
                            tmp=nirs.core.Data;
                            if obj.multitrial % 3-D stack of blocks to calc FC once
                                dur = min(s.dur) - 2*obj.ignore;
                                nsamp = floor( dur * data(idxA).Fs );
                                tmp.data = zeros([nsamp 2*size(dataA,2) length(s.onset)]);
                                for j=1:length(s.onset)
                                    if(obj.verbose)
                                        disp(['   ' num2str(j) ' of ' num2str(length(s.onset))]);
                                    end
                                    startpt = find(time>(s.onset(j)+obj.ignore),1,'first');
                                    lstpts=startpt:(startpt+nsamp-1);
                                    tmp.data(:,:,j)=[dataA(lstpts,:) dataB(lstpts,:)];
                                end
                                tmp.time=time(lstpts);
                                [r,p,dfe]=obj.corrfcn(tmp);
                                
                            else  % Iterate over each block to calc FC
                                for j=1:length(s.onset)
                                    if(obj.verbose)
                                        disp(['   ' num2str(j) ' of ' num2str(length(s.onset))]);
                                    end
                                    lstpts=find(time>s.onset(j)+obj.ignore &...
                                        time<s.onset(j)+s.dur(j)-obj.ignore);
                                    tmp.data=[dataA(lstpts,:) dataB(lstpts,:)];
                                    tmp.time=time(lstpts);
                                    [r(:,:,j),p,dfe(j)]=obj.corrfcn(tmp);
                                end
                            end
                            
                            if(obj.symetric)
                                r = atanh(r); % r-to-Z
                                if contains(func2str(obj.corrfcn),'nirs.sFC.grangers')
                                    r = exp(2*r); % Z-to-F
                                end
                                for j = 1:size(r,3)
                                    aa=r(1:end/2,1:end/2,j);            % within subject A
                                    ab=r(1:end/2,end/2+1:end,j);        % from A to B
                                    ba=r(end/2+1:end,1:end/2,j);        % from B to A
                                    bb=r(end/2+1:end,end/2+1:end,j);    % within subject B
                                    r(:,:,j) = ([aa ab; ba bb] + [bb ba; ab aa]) ./ 2;
                                end
                                if contains(func2str(obj.corrfcn),'nirs.sFC.grangers')
                                    r = log(r)/2; % F-to-Z
                                end
                                r = tanh(r); % Z-to-r
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
                    tmp=nirs.core.Data;
                    tmp.data=[dataA dataB];
                    tmp.time=time;
                    
                    lst=find(tmp.time<obj.ignore | tmp.time>tmp.time(end)-obj.ignore);
                    tmp.data(lst,:)=[];
                    tmp.time(lst)=[];
                    [r,p,dfe]=obj.corrfcn(tmp);
                    
                    if(obj.symetric)
                        r = atanh(r); % r-to-Z
                        if contains(func2str(obj.corrfcn),'nirs.sFC.grangers')
                            r = exp(2*r); % Z-to-F
                        end
                        aa=r(1:end/2,1:end/2);            % within subject A
                        ab=r(1:end/2,end/2+1:end);        % from A to B
                        ba=r(end/2+1:end,1:end/2);        % from B to A
                        bb=r(end/2+1:end,end/2+1:end);    % within subject B
                        r = ([aa ab; ba bb] + [bb ba; ab aa]) ./ 2;
                        if contains(func2str(obj.corrfcn),'nirs.sFC.grangers')
                            r = log(r)/2; % F-to-Z
                        end
                        r = tanh(r); % Z-to-r
                    end
                    
                    connStats(i).dfe=dfe;
                    connStats(i).R=r;
                    connStats(i).conditions=cellstr('rest');
                end
                
                
                
                % Compute data hash and load cached result if match is found
                %                 clear hash cache_file
                %                 if exist(obj.cache_dir,'dir')
                %                     hashopt.Method = 'SHA-256';
                %                     tmpobj = obj;
                %                     tmpobj.link = []; tmpobj.cache_dir = []; tmpobj.cache_rebuild = [];
                %                     tmpobj.corrfcn = func2str(tmpobj.corrfcn); tmpobj.prevJob = [];
                %                     hash = DataHash( {dataA,dataB,tmpobj,mask} , hashopt );
                %                     cache_file = fullfile( obj.cache_dir , [hash '.mat'] );
                %                     if ~obj.cache_rebuild && exist(cache_file,'file')
                %                         tmp = load(cache_file,'connStats');
                %                         connStats(i).R = tmp.connStats.R;
                %                         connStats(i).dfe = tmp.connStats.dfe;
                %                         disp(['Finished ' num2str(i) ' of ' num2str(height(obj.link)) ' (cached)'])
                %                         continue;
                %                     end
                %                 end
                
                connStats(i).R(1:end/2,1:end/2,:)=0;
                connStats(i).R(1+end/2:end,1+end/2:end,:)=0;
                
                disp(['Finished ' num2str(i) ' of ' num2str(height(obj.link))])
                
                %                 if exist('cache_file','var')
                %                     tmp = [];
                %                     tmp.connStats = connStats(i);
                %                     save(cache_file,'-struct','tmp');
                %                 end
                
            end
            
        end
    end
end

