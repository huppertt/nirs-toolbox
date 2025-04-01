classdef Hyperscanning_GLM < nirs.modules.AbstractGLM
    %% GLM- this is a wrapper for the OLS, AR-IRLS, and SPM-NIRS regression programs%
    %
    % Options:
    %     basis       - a Dictionary object containing temporal bases using stim name as key
    %     verbose     - flag to display progress
    %     trend_func  - a function that takes in a time vector and returns trend regressors
    %     type        - {OLS, NIRS-SPM, or [AR-IRLS]}
    % Example:
    %     j = nirs.modules.GLM();
    %     j.type = 'AR-IRLS';
    %     b = Dictionary();
    %     b('default') = nirs.design.basis.Canonical(); % default basis
    %     b('A')       = nirs.design.basis.Gamma();     % a different basis for condition 'A'
    %
    %     j.basis = b;
    %
    %     j.trend_func = @(t) nirs.design.trend.legendre(t, 3); % 3rd order polynomial for trend
    %
    % Note:
    %     trend_func must at least return a constant term unless all baseline periods are
    %     specified explicitly in the stimulus design with BoxCar basis functions
    properties
        type;
        method = 'Pearson';  % WaveletMag([.05 .1]), WaveletPhase([.05 .1]);
        AddShortSepRegressors = false;
        link;
        linkVariable = 'hyperscan'; % hyperscan variable from nirx
        options;
    end
    methods
        function obj = Hyperscanning_GLM( prevJob )
            if nargin > 0, obj.prevJob = prevJob; end
            
            obj.name = 'GLM model';
            obj.basis('default') = nirs.design.basis.Canonical();
            obj.type='AR-IRLS';
            obj.options=[];
        end
        
        function obj = set.type(obj,type)
            validtypes={'OLS','NIRS-SPM','AR-IRLS','MV-GLM','Nonlinear','Ordinal','RepeatedMeasures'};
            if(~ismember(type,validtypes))
                disp('type must be one of : ')
                disp(strvcat(validtypes));
                return;
            else
                obj.type=type;
            end
            % use the call functions to evoke any special messages or
            % conditions (e.g. NIRS-SPM checks for code in path);
            switch(obj.type)
                case('OLS')
                    j=nirs.modules.OLS();
                case('AR-IRLS');
                    j=nirs.modules.AR_IRLS();
                case('NIRS-SPM')
                    j=nirs.modules.NIRS_SPM_GLM();
                case('MV-GLM')
                    j=nirs.modules.MultiVarGLM();
                    disp(['Note: Inputs expected to be optical density']);
                case('Nonlinear')
                    j=nirs.modules.nonlin_GLM();
                    obj.basis=j.basis;
                case('Ordinal')
                    j=nirs.modules.MMR_GLM();
                case('RepeatedMeasures')
                    j=nirs.modules.repeatedMeas_GLM();
                otherwise
                    error('type not recognized');
            end
            obj.citation=j.citation;
            obj.options=[];
            pj=obj.prevJob;
            flds=fields(j);
            lst=find(ismember(flds,fields(obj)));
            for i=1:length(lst)
                obj.(flds{lst(i)})=j.(flds{lst(i)});
            end
            lst=find(~ismember(flds,fields(obj)));
            for i=1:length(lst)
                obj.options=setfield(obj.options,flds{lst(i)},j.(flds{lst(i)}));
            end
            obj.prevJob=pj;
        end
        
        function S = runThis( obj, data )
             
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

                    for i=1:length(uniquefiles)
                        lst=find(ib==i);
                        ScanA(i,1)=lst(1);
                        ScanB(i,1)=lst(2);
                        if(length(obj.linkVariable)>1)
                            relationship{i,1}=tbl(lst,:).(obj.linkVariable{2});
                        else
                            relationship{i,1}=[];
                        end
                    end

                    OffsetA = zeros(size(ScanA));  % The time shift of the "A" files (in sec)
                    OffsetB = zeros(size(ScanB));  % The time shift of the "B" files (in sec)

                    link = table(ScanA,ScanB,OffsetA,OffsetB,relationship);
                    obj.link=link;
                else
                    warning('link variable must be specified');
                    S=nirs.core.ChannelStats;
                    return
                end
            end

            if (obj.AddShortSepRegressors)
                j=nirs.modules.AddShortSeperationRegressors();
                j=nirs.modules.RemoveShortSeperations(j);
            else
                j=nirs.modules.Assert;
                j.condition=@(data)isa(data,'nirs.core.Data');
            end
            data=j.run(data);
            
            switch(obj.type)
                case('OLS')
                    j=nirs.modules.OLS;
                case('AR-IRLS');
                    j=nirs.modules.AR_IRLS;
                case('NIRS-SPM')
                    j=nirs.modules.NIRS_SPM_GLM;
                case('MV-GLM')
                    j=nirs.modules.MultiVarGLM;
                case('Nonlinear')
                    j=nirs.modules.nonlin_GLM;
                case('Ordinal')
                    j=nirs.modules.MMR_GLM;
                case('RepeatedMeasures')
                    j=nirs.modules.repeatedMeas_GLM;
                otherwise
                    error('type not recognized');
            end
            j.basis=obj.basis;
            j.verbose=obj.verbose;
            j.trend_func=obj.trend_func;
            j.goforit=obj.goforit;
            
            
            if(~isempty(obj.options))
                flds=fields(obj.options);
                for i=1:length(flds)
                    j.(flds{i})=obj.options.(flds{i});
                end
            end
            
            
            if (obj.AddShortSepRegressors)
                Stim=unique(nirs.getStimNames(data));
                for ii=1:length(Stim)
                    Stim{ii}=[Stim{ii} '*'];
                end
                
                j=nirs.modules.KeepStims(j);
                j.regex=true;
                j.listOfStims=Stim;
            end

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
                    time=[max(timeA(1),timeB(1)):1/data(idxA).Fs:min(timeA(end),timeB(end))]';
                    for id=1:size(dataA,2)
                        dataA(1:length(time),id)=interp1(timeA,dataA(:,id),time);
                    end
                    for id=1:size(dataB,2)
                        dataB(1:length(time),id)=interp1(timeB,dataB(:,id),time);
                    end
                    dataA=dataA(1:length(time),:);
                    dataB=dataB(1:length(time),:);
                end

                data2(i)=make_ts(dataA,dataB, ...
                    data(idxA).probe,data(idxB).probe,data(idxA).Fs,obj.method);
                data2(i).time=data(idxA).time;
                data2(i).stimulus=data(idxA).stimulus;
                keys=data(idxB).stimulus.keys;
                for ii=1:length(keys)
                    if(contains(keys{ii},'SS_PCA'))
                        data2(i).stimulus(strrep(keys{ii},'PCA','PCA_B'))=data(idxB).stimulus(keys{ii});
                    end
                end

                data2(i).demographics=data(idxA).demographics; % For now



            end

                S=j.run(data2);
            
            for idx=1:length(S)
                [~,Stim]=nirs.design.createDesignMatrix(data(idx).stimulus,data(idx).time,obj.basis);
                StimNew=unique(nirs.getStimNames(S(idx)));
                
                lst=[]; lst2=[];
                for j=1:length(Stim)
                    ss=Stim{j};
                    if(~isempty(strfind(ss,':')))
                        ss(strfind(ss,':'):end)=[];
                    end
                    
                    st=data(idx).stimulus(ss);
                    if(~ismember('regressor_no_interest',fields(st)))
                        lst=[lst j];
                    else
                        if(~st.regressor_no_interest)
                            
                            lst=[lst j];
                        else
                            lst2=[lst2 j];
                        end
                    end
                end
                StimRm={Stim{lst2}};
                Stim={Stim{lst}};
                
                SS={};
                for i=1:length(StimNew)
                    for j=1:length(Stim)
                        if(~isempty(ismember(StimNew{i},Stim{j},'rows')) & ...
                                isempty(strfind(StimNew{i},'SS_PCA')) & ...
                                ~ismember(StimNew{i},StimRm))
                            SS{end+1}=StimNew{i};
                        end
                    end
                end
                SS=unique(SS);
                for ii=1:length(SS)
                    SS{ii}=[SS{ii} '*'];
                end
                
                j=nirs.modules.KeepStims;
                j.listOfStims=SS;
                j.regex=true;
                S(idx)=j.run(S(idx));
            end
            
            
            
            
            
            
        end
        
    end
    
end

function data2= make_ts(dataA,dataB,probeA,probeB,Fs,method)

cnt=1;
link=struct;
if(strcmp(method,'Pearson'))
    for a=1:size(dataA,2)
        for b=1:size(dataB,2)
            dd(:,cnt)=dataA(:,a).*dataB(:,b);
            cnt=cnt+1;
        end
    end
elseif(contains(method,'Wavelet'))
    str=strrep(method,'Wavelet','');
    str=strrep(str,'(','');
    str=strrep(str,')','');
    str=strrep(str,'Mag','');
    str=strrep(str,'Phase','');
    freqband = str2num(str);
    freqband=freqband/Fs;
    for i=1:size(dataA,2);
        [wtA(:,:,i),f]=cwt(dataA(:,i));
    end;
    for i=1:size(dataB,2);
        [wtB(:,:,i),f]=cwt(dataB(:,i));
    end;

    cnt=1;
    for i=1:size(dataA,2);
        for j=1:size(dataB,2);
            WCOH(:,:,cnt)=wtA(:,:,i).*conj(wtB(:,:,j));
            cnt=cnt+1;
        end
    end
    

    dd = zeros(size(WCOH,2),size(WCOH,3));
   
    lst=find(f>=freqband(1) & f<=freqband(2));
    if(contains(method,'Mag'))
        dd = squeeze(median(abs(WCOH(lst,:,:)),1));
    else
        dd= squeeze(median(angle(WCOH(lst,:,:)),1));
    end
    
end

cnt=1;
lst=[];
for a=1:size(dataA,2)
    for b=1:size(dataB,2)
        link.source{cnt,1}=['src-' num2str(probeA.link.source(a)) ':' ...
            'det-' num2str(probeA.link.detector(a))];
        link.detector{cnt,1}=['src-' num2str(probeB.link.source(b)) ':' ...
            'det-' num2str(probeB.link.detector(b))];
        if(iscell(probeA.link.type))
            if(probeA.link.type{a}==probeB.link.type{b})
                link.type{cnt,1}=probeA.link.type{a};
               
            else
                lst=[lst; cnt];
                link.type{cnt,1}=' ';
            end
        else
            if(probeA.link.type(a)==probeB.link.type(b))
                link.type(cnt,1)=probeA.link.type(a);
            else
                lst=[lst; cnt];
                link.type(cnt,1)=NaN;
            end

        end
        cnt=cnt+1;
    end;
end;
data2=nirs.core.Data;
data2.probe=probeA;
data2.data=dd;
data2.probe.link=struct2table(link);
data2.probe.link(lst,:)=[];
data2.data(:,lst)=[];

end