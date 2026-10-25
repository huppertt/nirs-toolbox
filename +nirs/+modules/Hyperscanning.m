classdef Hyperscanning < nirs.modules.AbstractModule
    %% Hyper-scanning - Computes all-to-all connectivity model between two seperate files
    % Outputs nirs.core.ConnStats object
    
    properties
        corrfcn;  % function to use to compute correlation (see +nirs/+sFC for options)
        divide_events;  % if true will parse into multiple conditions
        min_event_duration;  % minimum duration of events
        link;
        linkVariable = 'hyperscan'; % hyperscan variable from nirx
        symmetric;
        verbose;
        ignore;  % seconds at start/end of each scan or block to ignore
        estimate_null; % flag to also estimate connectivity between all non-paired subjects 
                       %  and creates an 'Pairing' demographic field in the output with values 'Actual' or 'Null'
    end

    methods
        function obj = Hyperscanning( prevJob )
            obj.name = 'Hyperscanning';
            obj.corrfcn = @(data)nirs.sFC.ar_corr(data,'18xFs',true,10);  %default to use AR-whitened robust correlation (BIC capped at 10)
            obj.divide_events=false;
            obj.min_event_duration=30;
            obj.symmetric=true;
            obj.ignore=10;
            obj.verbose=false;
            obj.estimate_null=false;
            
            obj.link=table(1,1,0,0,'VariableNames',{'ScanA','ScanB','OffsetA','OffsetB'});
            obj.link(1,:)=[];
            
            if nargin > 0
                obj.prevJob = prevJob;
            end
            
            obj.citation='Santosa, H., Aarabi, A., Perlman, S. B., & Huppert, T. J. (2017). Characterization and correction of the false-discovery rates in resting state connectivity using functional near-infrared spectroscopy. Journal of Biomedical Optics, 22(5), 055002-055002.';
            
        end
        
        function connStats = runThis( obj, data )
            if(~iscell(obj.linkVariable))
                obj.linkVariable={obj.linkVariable};
            end
            
            

            if(isempty(obj.link))
                % NIRx files have a hyperscan variable upon loading that I
                % can use here
                tbl=nirs.createDemographicsTable(data);
                [tbl,idx]=sortrows(tbl,obj.linkVariable);
                data=data(idx);
                

                if(ismember(obj.linkVariable{1},nirs.createDemographicsTable(data).Properties.VariableNames))
                    hyperscanfiles=nirs.createDemographicsTable(data).(obj.linkVariable{1});
                    
                    for i=1:length(hyperscanfiles); if(isempty(hyperscanfiles{i})); hyperscanfiles{i}=''; end; end;
                    uniquefiles=unique(hyperscanfiles);
                    [ia,ib]=ismember(hyperscanfiles,uniquefiles);
                    
                    cnt=1;
                    for i=1:length(uniquefiles)
                        lst=find(ib==i);
                        if(length(lst)==2)
                            ScanA(cnt,1)=lst(1);
                            ScanB(cnt,1)=lst(2);
                            if(length(obj.linkVariable)>1)
                                relationship{cnt,1}=tbl(lst,:).(obj.linkVariable{2});
                            else
                                relationship{cnt,1}=[];
                            end
                            cnt=cnt+1;
                        end
                    end
                    
                    OffsetA = zeros(size(ScanA));  % The time shift of the "A" files (in sec)
                    OffsetB = zeros(size(ScanB));  % The time shift of the "B" files (in sec)
                    
                    link = table(ScanA,ScanB,OffsetA,OffsetB,relationship);
                    obj.link=link;
                else
                    warning('link variable must be specified');
                    connStats=nirs.core.sFCStats();
                    return
                end
            end
            
            % Add null distribution pairings
            if obj.estimate_null
                
                obj.link.isNull = false(height(obj.link),1);
                
                % Create list of all file pairings
                combos = combnk(1:length(data),2);
                
                % Remove true pairings from the list
                combos = setdiff(combos,[obj.link.ScanA obj.link.ScanB],'rows');
                
                % Check the number is correct
                num_nulls = (length(data)^2 - 2*length(data))/2;
                assert(size(combos,1)==num_nulls,'Unexpected number of null pairings: %i',size(combos,1));
                
                % Create null link table
                link_null = table(combos(:,1),combos(:,2),zeros(num_nulls,1),zeros(num_nulls,1),true(num_nulls,1),...
                    'VariableNames',{'ScanA','ScanB','OffsetA','OffsetB','isNull'});
                
                obj.link = [obj.link; link_null];
                
            end
            
            nDyads = height(obj.link);
            connStats = repmat(nirs.core.sFCStats(), 1, nDyads);

            % Cache loop-invariant values for parfor
            hasRelationship = ismember('relationship', obj.link.Properties.VariableNames);
            linkTbl = obj.link;
            corrfcn_   = obj.corrfcn;
            divide_    = obj.divide_events;
            symmetric_ = obj.symmetric;
            ignore_    = obj.ignore;
            minDur_    = obj.min_event_duration;
            estNull_   = obj.estimate_null;
            verbose_   = obj.verbose;

            parfor i = 1:nDyads

                idxA = linkTbl.ScanA(i);
                idxB = linkTbl.ScanB(i);

                dataA = data(idxA).data;
                timeA = data(idxA).time + linkTbl.OffsetA(i);
                dataB = data(idxB).data;
                timeB = data(idxB).time + linkTbl.OffsetB(i);

                % Make sure we are using the same time base
                if isequal(timeA, timeB)
                    time = timeA;
                else
                    time = (max(timeA(1),timeB(1)):1/data(idxA).Fs:min(timeA(end),timeB(end)))';
                    for id = 1:size(dataA,2)
                        dataA(1:length(time),id) = interp1(timeA, dataA(:,id), time);
                    end
                    dataB_orig = dataB;
                    for id = 1:size(dataB,2)
                        dataB(1:length(time),id) = interp1(timeB, dataB_orig(:,id), time);
                    end
                    dataA = dataA(1:length(time),:);
                    dataB = dataB(1:length(time),:);
                end

                cs = nirs.core.sFCStats();
                cs.type = corrfcn_;
                cs.description = data(idxA).description;

                if hasRelationship
                    relationship = linkTbl(i,:).relationship{1};
                    if isa(data(idxA).probe, 'nirs.core.Probe1020')
                        cs.probe = nirs.core.ProbeHyperscan1020([data(idxA).probe, data(idxB).probe], relationship);
                    else
                        cs.probe = nirs.core.ProbeHyperscan([data(idxA).probe, data(idxB).probe], relationship);
                    end
                else
                    if isa(data(idxA).probe, 'nirs.core.Probe1020')
                        cs.probe = nirs.core.ProbeHyperscan1020([data(idxA).probe, data(idxB).probe]);
                    else
                        cs.probe = nirs.core.ProbeHyperscan([data(idxA).probe, data(idxB).probe]);
                    end
                end

                demog = {data(idxA).demographics; data(idxB).demographics};
                if estNull_
                    if linkTbl.isNull(i)
                        demog{1}('Pairing') = {'Null'};
                        demog{2}('Pairing') = {'Null'};
                    else
                        demog{1}('Pairing') = {'Actual'};
                        demog{2}('Pairing') = {'Actual'};
                    end
                end
                cs.demographics = demog;
                cs.R = [];

                if divide_
                    stim = data(idxA).stimulus;

                    cnt = 1;
                    for idx = 1:length(stim.keys)
                        s = stim(stim.keys{idx});
                        lst = find(s.dur - 2*ignore_ > minDur_);
                        if ~isempty(lst)
                            s.onset = s.onset(lst);
                            s.dur = s.dur(lst);

                            n1 = size(dataA,2) + size(dataB,2);
                            r = zeros(n1, n1, length(s.onset));
                            dfe = zeros(length(s.onset), 1);

                            for j = 1:length(s.onset)
                                tmp = nirs.core.Data;
                                lstpts = find(time > s.onset(j)+ignore_ & ...
                                    time < s.onset(j)+s.dur(j)-ignore_);
                                tmp.data = [dataA(lstpts,:) dataB(lstpts,:)];
                                tmp.time = time(lstpts);
                                [r(:,:,j), ~, dfe(j)] = corrfcn_(tmp);
                            end

                            if symmetric_
                                r = atanh(r);
                                for j = 1:size(r,3)
                                    aa = r(1:end/2, 1:end/2, j);
                                    ab = r(1:end/2, end/2+1:end, j);
                                    ba = r(end/2+1:end, 1:end/2, j);
                                    bb = r(end/2+1:end, end/2+1:end, j);
                                    r(:,:,j) = ([aa ab; ba bb] + [bb ba; ab aa]) ./ 2;
                                end
                                r = tanh(r);
                            end

                            cs.dfe(cnt) = sum(dfe);
                            cs.R(:,:,cnt) = tanh(mean(atanh(r), 3));
                            cs.conditions{cnt} = stim.keys{idx};
                            cnt = cnt + 1;
                        end
                    end

                else
                    tmp = nirs.core.Data;
                    tmp.data = [dataA dataB];
                    tmp.time = time;

                    lst = find(tmp.time < ignore_ | tmp.time > tmp.time(end) - ignore_);
                    tmp.data(lst,:) = [];
                    tmp.time(lst) = [];
                    [r, ~, dfe] = corrfcn_(tmp);

                    if symmetric_
                        r = atanh(r);
                        if ~isempty(strfind(func2str(corrfcn_), 'nirs.sFC.grangers'))
                            r = exp(2*r);
                        end
                        aa = r(1:end/2, 1:end/2);
                        ab = r(1:end/2, end/2+1:end);
                        ba = r(end/2+1:end, 1:end/2);
                        bb = r(end/2+1:end, end/2+1:end);
                        r = ([aa ab; ba bb] + [bb ba; ab aa]) ./ 2;
                        if ~isempty(strfind(func2str(corrfcn_), 'nirs.sFC.grangers'))
                            r = log(r)/2;
                        end
                        r = tanh(r);
                    end

                    cs.dfe = dfe;
                    cs.R = r;
                    cs.conditions = cellstr('rest');
                end

                connStats(i) = cs;

            end
            
        end
    end
end

