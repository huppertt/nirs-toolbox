classdef repeatedMeas_GLM < nirs.modules.AbstractGLM
%% AR_IRLS - Performs first-level per file GLM analysis.
% 
% Options:
%     basis       - a Dictionary object containing temporal bases using stim name as key
%     verbose     - flag to display progress
%     trend_func  - a function that takes in a time vector and returns trend regressors
%     
% Example:
%     j = nirs.modules.AR_IRLS();
%     
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
        formula = 'Time';
     end
    
     
    methods
        function obj = repeatedMeas_GLM( prevJob )
            if nargin > 0, obj.prevJob = prevJob; end
            
            obj.name = 'GLM with repeated measurements model';
            obj.basis('default') = nirs.design.basis.Canonical();
%             obj.citation=['Barker, Jeffrey W., Ardalan Aarabi, and Theodore J. Huppert.'...
%                 '"Autoregressive model based algorithm for correcting motion and serially '...
%                 'correlated errors in fNIRS." Biomedical optics express 4.8 (2013): 1366-1379.'];
            
        end
        
        function S = runThis( obj, data )
            vec = @(x) x(:);
            
            for i = 1:numel(data)
                % get data
                d  = data(i).data;
                t  = data(i).time;
                Fs = data(i).Fs;
                
                probe = data(i).probe;
                
                % make sure data is in order
                
                if(~isempty(strfind(class(probe),'nirs')))
                    if(~ismember('source',probe.link.Properties.VariableNames) & ...
                    ismember('ROI',probe.link.Properties.VariableNames))
                        [probe.link, idx] = nirs.util.sortrows(probe.link, {'ROI','type'});
                    else
                        [probe.link, idx] = nirs.util.sortrows(probe.link, {'source', 'detector','type'});
                    end
                elseif(~isempty(strfind(class(probe),'eeg')))
                    [probe.link, idx] = nirs.util.sortrows(probe.link, {'electrode','type'});
                else
                    error('data type not supported');
                end
                    d = d(:, idx);
                    % get experiment design
                    C = obj.getTrendMatrix( t );
                    
                    a=obj.basis('default');
                    a1 = a.peakTime;
                    a2 = a.uShootTime;
                    b1 = a.peakDisp;
                    b2 = a.uShootDisp;
                    c  = a.ratio;
                    
                    % sampling freq
                   t2 = (0:1/Fs:a.duration)';
                    d(end+[1:length(t2)],:)=NaN;
                    IRF = a.getImpulseResponse( a1, b1, a2, b2, c, t2 );
                    
                    for chIdx=1:size(d,2)
                        
                        Y=zeros(length(t2),1);
                        cnt=1;
                        
                        
                        cond={};
                        keys=data(i).stimulus.keys;
                        metadata=[];
                        for ii=1:length(keys)
                            st=data(i).stimulus(keys{ii});
                            for j=1:length(st.onset)
                                [~,id]=min(abs(t-st.onset(j)));
                                Y(:,cnt)=d(id+[0:length(t2)-1],chIdx);
                                cond{cnt,1}=keys{ii};
                                cnt=cnt+1;
                            end
                            metadata=[metadata; st.metadata];
                        end
                        within=[table(cond) metadata];
                        within.Time=[1:height(within)]';
                        
                        between=struct;
                        for ii=1:length(cond)
                            between=setfield(between,['y' num2str(ii)],Y(:,ii));
                        end
                        str=['y1-y' num2str(length(cond)) '~1+'];
                       str2=['y1-y' num2str(length(cond)) '~1'];
                       
                        for jj=1:length(keys)
                        for ii=1:size(IRF,2)
                           if(size(IRF,2)>1)
                                between=setfield(between,[keys{jj} '_' num2str(ii)],IRF(:,ii));
                                str=[str  keys{jj} '_' num2str(ii) '+'];
                           else
                                between=setfield(between,[keys{jj} ],IRF(:,ii));
                                str=[str keys{jj} '+' ];
                           end
                        end
                        end
                        str(end)=[];
                        r{chIdx}=nirs.math.fitrm(struct2table(between),str,1,fix(data(i).Fs*4),'WithinDesign',within);
                        
                        fprintf(1,'%d of %d\n',chIdx,size(d,2))
                    end
                     
                
                    var=[];

                    S(i)=nirs.core.ChannelFStats;
                    S(i).description=['Repeated ANOVA: ' obj.formula];
                    S(i).probe=data(i).probe;
                    S(i).demographics=data(i).demographics;
                    
                    
                    
                    for chIdx=1:length(r)
                        tbl=ranova(r{chIdx},'WithinModel',obj.formula);
                        for ii=1:length(keys)
                            for jj=1:height(tbl)
                                if(~isempty(strfind(tbl.Row{jj},[keys{ii}])))
                                    n=tbl(jj,:);
                                    nn=n.Properties.RowNames;
                                    n.Properties.RowNames={};
                                    
                                    var=[var; [table(nn,'VariableNames',{'cond'}) S(i).probe.link(chIdx,:) n table({r{chIdx}},'VariableNames',{'model'})]];
                                end
                            end
                        end
                        a=anova(r{chIdx});
                        
                          
                    end
                    S(i).F=var.F;
                    S(i).df2=var.DF(1);
                    S(i).df1=length(d)-length(keys)-1;
                    
                    var.F=[];
                    var.DF=[];
                    S(i).variables=var;
                
                % print progress
                if(obj.verbose)
                 obj.printProgress( i, length(data) )
                end
            end

        end
        
        
        function prop = javaoptions(obj)
            
            prop=javaoptions@nirs.modules.AbstractGLM(obj);
            opts=obj.options;
            
            diction=nirs.util.createDictionaryFromToolset('nirs.design.basis');
            DictionaryProp=javatypes('enum',{diction.values});
            set(DictionaryProp,'Name','basis','Value','test');
            set(DictionaryProp,'Category','Misc');
            set(DictionaryProp,'Description','Select the canonical basic function');
            prop(find(ismember(opts,'basis')))=DictionaryProp;
  
            
        end
        
    end
    
end

