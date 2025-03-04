classdef ChannelStats
    %% CHANNELSTATS - Holds regression stats for channel space.
    %
    % Properties:
    %     description  - description of data (e.g. filename)
    %     variables    - a table containing w/ source, detector, type, cond for
    %                    each regression coefficient
    %     beta         - regression coefficients
    %     covb         - covariance of beta
    %     dfe          - degrees of freedom
    %     probe        - Probe object describing measurement geometry
    %     demographics - Dictionary containing demographics info
    %                    (e.g. demographics('age') returns 28)
    %     conditions   - (dependent) list of stim conditions
    %     tstat        - (dependent) t-stats of beta
    %     p            - (dependent) p-values of beta
    %     q            - (dependent) q-values of beta (false discovery rate)
    %
    %  Methods:
    %     getCritT    - returns critical t value
    %     draw        - draws beta or tstat values on top of probe geometry
    %     table       - returns a table of all stats (minus full covariance)
    %     ttest       - performs a t-test and returns a new ChannelStats object
    %     ftest       - performs an F-test across conditions
    %     jointTest   - performs a joint hypothesis test across all channels in
    %                   each SD pair (i.e. joint test of hbo & hbr for S1-D1)
    %     sorted      - returns new ChannelStats object with sorted channels
    
    properties
        description     % description of data (e.g. filename)
        
        variables       % table describing regression coefficients
        
        beta            % regression coefficients
        covb            % covariance of beta
        dfe          	% degrees of freedom
        tail ='two-sided'; % {two-sided,right,left}  tail for ttest
        probe           % Probe object describing measurement geometry
        demographics    % Dictionary containing demographics info
        basis           % basis set info used to create model
    end
    
    properties ( Dependent = true )
        conditions      % (dependent) list of stim conditions
        
        tstat           % (dependent) t-stats of beta
        p               % (dependent) p-values of beta
        
        q               % (dependent) q-values of beta (false discovery rate)
    end
    properties (Hidden = true )
        pvalue_fixed;
        categoricalvariableInfo;
        WhiteningW;
        tag;
        UserData;
    end
    
    
    
    
    
    
    methods
        
        
        
        
        function obj = set.tail(obj,tail)
            validtypes={'left','right','one-sided','two-sided'};
            if(~ismember(lower(tail),validtypes))
                disp('type must be one of : ')
                disp(strvcat(validtypes));
                return;
            else
                obj.tail=lower(tail);
            end
        end
        
        function hrf=HRF(obj,type)
            % function extracts the HRF from the stats variable
            
            if(nargin<2)
                type='hrf';
            end
            hrf=[];
            for i=1:length(obj)
                if(isempty(hrf))
                    hrf=nirs.design.extractHRF(obj(i),obj(i).basis.base,obj(i).basis.stim,obj(i).basis.Fs,type);
                else
                    hrf=[hrf;nirs.design.extractHRF(obj(i),obj(i).basis.base,obj(i).basis.stim,obj(i).basis.Fs,type)];
                end
            end
        end
        
        % unique conditions
        function c = get.conditions( obj )
            if ~isempty(obj.variables)
                c = unique(obj.variables.cond);
            else
                c = [];
            end
        end
        
        % t statistic calculation
        function tstat = get.tstat( obj )
            tstat = obj.beta ./ sqrt(diag(obj.covb));
        end
        
        function p = get.p( obj )
            
            if(~isempty(obj.pvalue_fixed))
                p=obj.pvalue_fixed;  % hack to allow reuse in the Conn Seed model
                return
            end
            
            t = obj.tstat;
            
            if(strcmp(lower(obj.tail),'right') | strcmp(lower(obj.tail),'one-sided'))
                p = tcdf(-t, obj.dfe);
            elseif(strcmp(lower(obj.tail),'left'))
                p = tcdf(t, obj.dfe);
            else % two-side or unknown
                p = 2*tcdf(-abs(t), obj.dfe);
            end
        end
        
        % q values
        function q = get.q( obj )
            q = reshape( nirs.math.fdr( obj.p(:) )', size(obj.p) );
        end
        
        
        function scatterplot(obj,channelIdx,confinterval,seperateplot)
            
            if(nargin<3 || isempty(confinterval))
                confinterval=.95;
            end
            if(nargin<4)
                seperateplot=false;
            end
            
            
            
            
            var = obj.variables;
            if(~ismember('model',var.Properties.VariableNames))
                disp('WARNING: The Stats variable does not contain the linear model info')
                disp('Either this is a 1st-level statistic or the group model was run');
                disp('without the "include_diagnostics" flag');
                return
            end
            if(nargin<2)
                warning('Channel Index [which row in the probe.link variable] is required');
                return
            end
            ch=obj.probe.link(channelIdx,:);
            
            
            if(isa(obj.probe,'nirs.core.ProbeROI'))
                lst=find(ismember(var.ROI,ch.ROI) &...
                ismember(var.type,ch.type));
            else
            lst=find(ismember(var.source,ch.source) & ismember(var.detector,ch.detector) &...
                ismember(var.type,ch.type));
            end
            
            % all models for that src/det/type will be the same
            model = var.model{lst(1)};
            modelvals=model.Variables;
            
            
            m=abs(tinv((1-confinterval)/2,model.DFE));
            
            if(~isempty(obj.WhiteningW))
                iW=inv(obj.WhiteningW);
                flds=modelvals.Properties.VariableNames;
                mask=zeros(size(iW,1),1);
                for i=1:length(flds)
                    modelvals.(flds{i})=full(iW*modelvals.(flds{i}));
                    
                    if(~strcmp(flds{i},model.ResponseName))
                        mask=max(mask,abs(modelvals.(flds{i})));
                    end
                end
                EB=abs(sqrt(diag(obj.WhiteningW)));
                modelvals(find(abs(mask)<1E-10),:)=[];
                EB(find(abs(mask)<1E-10))=[];
            else
                EB=[];
            end
            
            categoricalvariableInfo=obj.categoricalvariableInfo;
            variablename=categoricalvariableInfo.Properties.RowNames;
            
            
            nNonCat = length(find(~categoricalvariableInfo.IsCategorical));
            nNonCat=0;
             lstNonCat={};
            for i=1:length(obj.conditions)
                 s=strsplit(obj.conditions{i},':');
                if(any(ismember(s,{categoricalvariableInfo.Row{~categoricalvariableInfo.IsCategorical}})))
                    nNonCat=nNonCat+1;
                     lstNonCat{end+1,1}=obj.conditions{i};
                end
            end
            
            
            
            
            m2=nNonCat+1;
            mm=1;
            
            % first plot the catgorical
            
            lst=find(categoricalvariableInfo.IsCategorical);
            for id=1: height(categoricalvariableInfo)
                if(~categoricalvariableInfo.IsCategorical(id))
                    if(categoricalvariableInfo.Range{id}(1)==0 & categoricalvariableInfo.Range{id}(2)==1)
                        lst=[lst id];
                    end
                end
            end

            if(~isempty(lst))
                if(categoricalvariableInfo.IsCategorical(lst(1)))

                    name=categoricalvariableInfo.Range{lst(1)}';
                    if(~strcmp(categoricalvariableInfo.Row{1},'cond'))
                        for ii=1:length(name)
                            name{ii}=[categoricalvariableInfo.Row{1} '_' name{ii}];
                        end
                    end
                else
                    name=cellstr(categoricalvariableInfo.Row{lst(1)});
                end
                for i=2:length(lst)
                    if(categoricalvariableInfo.IsCategorical(lst(i)))
                        n2=categoricalvariableInfo.Range{lst(i)}';
                        if(~strcmp(categoricalvariableInfo.Row{i},'cond'))
                            for ii=1:length(n2)
                                n2{ii}=[categoricalvariableInfo.Row{i} '_' n2{ii}];
                            end
                        end
                    
                    name2=repmat(n2,length(name),1);


                    name=repmat(name,length(n2),1);
                    else
                        name2=cellstr(categoricalvariableInfo.Row{lst(i)});
                    end
                    name=horzcat(name,name2);

                end
            end

            if(~seperateplot)
                figure;
                set(groot,'defaultTextInterpreter','none')
                subplot(mm,m2,1);
                hold on;
                a=-999; b=999;
            else
                figure;
                set(groot,'defaultTextInterpreter','none')
                hold on;
                a=-999; b=999;
            end

            cnt=1; xlabels={};
            for i=1:length(obj.conditions)
                
                s=strsplit(obj.conditions{i},':');
                
                found=false;
                for j=1:size(name,1)
                    if(all(ismember(s,{name{j,:}})))
                        found=true;
                        
                    end
                end
                if(~found)
                   
                    continue;
                end
                
                beta=modelvals.beta;
                lstN=find(~ismember(model.Coefficients.Row,matlab.lang.makeValidName(obj.conditions{i})));
                for j=1:length(lstN)
                    beta=beta-model.Coefficients.Estimate(lstN(j))*modelvals.(model.Coefficients.Row{lstN(j)});
                end
                
                
                lst=find(ismember(model.Coefficients.Row,matlab.lang.makeValidName(obj.conditions{i})));
                xlabels{cnt}=obj.conditions{i};
                lst2=find(abs(modelvals.(model.Coefficients.Row{lst}))>eps(1));
                
                if(~isempty(EB))
                    e=nirs.util.scattererrorbar(cnt-1+1*(abs(modelvals.(model.Coefficients.Row{lst})(lst2))>1E-9),beta(lst2),EB(lst2));
                    set(e,'linewidth',2,'color','k');
                end
                s=scatter(cnt-1+1*(abs(modelvals.(model.Coefficients.Row{lst})(lst2))>1E-9),beta(lst2),'filled');
                l2(1)=line(cnt-1+[.75 1.25],...
                    [model.Coefficients.Estimate(lst)-model.Coefficients.SE(lst)...
                    model.Coefficients.Estimate(lst)-model.Coefficients.SE(lst)],'color','k');
                l2(2)=line(cnt-1+[.75 1.25],...
                    [model.Coefficients.Estimate(lst)+model.Coefficients.SE(lst)...
                    model.Coefficients.Estimate(lst)+model.Coefficients.SE(lst)],'color','k');
                l2(3)=line(cnt-1+[.75 .75],...
                    [model.Coefficients.Estimate(lst)+model.Coefficients.SE(lst)...
                    model.Coefficients.Estimate(lst)-model.Coefficients.SE(lst)],'color','k');
                l2(4)=line(cnt-1+[1.25 1.25],...
                    [model.Coefficients.Estimate(lst)+model.Coefficients.SE(lst)...
                    model.Coefficients.Estimate(lst)-model.Coefficients.SE(lst)],'color','k');
                
                l2(5)=line(cnt-1+[1 1],...
                    [model.Coefficients.Estimate(lst)+m*model.Coefficients.SE(lst)...
                    model.Coefficients.Estimate(lst)-m*model.Coefficients.SE(lst)],'color','k');
                l2(6)=line(cnt-1+[.85 1.15],...
                    [model.Coefficients.Estimate(lst)-m*model.Coefficients.SE(lst)...
                    model.Coefficients.Estimate(lst)-m*model.Coefficients.SE(lst)],'color','k');
                l2(7)=line(cnt-1+[.85 1.15],...
                    [model.Coefficients.Estimate(lst)+m*model.Coefficients.SE(lst)...
                    model.Coefficients.Estimate(lst)+m*model.Coefficients.SE(lst)],'color','k');
                l=line(cnt-1+[.75 1.25],...
                    [model.Coefficients.Estimate(lst) model.Coefficients.Estimate(lst)]);
                
                a=max(a,model.Coefficients.Estimate(lst)+m*model.Coefficients.SE(lst));
                b=min(b,model.Coefficients.Estimate(lst)-m*model.Coefficients.SE(lst));
                a=max(a,max(beta));
                b=min(b,min(beta));
                
                set(l,'LineWidth',4,'Color','r','LineStyle','-')
                set(l2,'LineWidth',2,'Color','k','LineStyle','-')
                
                set(s,'SizeData',150,'MarkerFaceColor','k','MarkerEdgeColor','k')
                cnt=cnt+1;
                
                
            end
            set(gca,'xlim',[.5 cnt-.5]);
            r=(a-b)*.05;
            
            set(gca,'ylim',[ b-r a+r]);
            ylabel(model.ResponseName)
            set(gca,'XTick',[1:cnt-1])
            set(gca,'XTickLabelRotation',45)
            for i=1:length(xlabels); xlabels{i}(find(xlabels{i}(:)=='_'))=[]; end
            set(gca,'XTickLabel',xlabels)
            set(gca,'fontsize',14)
            
            for i=1:length(lstNonCat)
                beta=modelvals.beta;
                lstN=find(~ismember(model.Coefficients.Row,lstNonCat{i}));
                LE=0;
                UE=0;
                 name=matlab.lang.makeValidName(lstNonCat{i});
                x=[min(modelvals.(name))-10:max(modelvals.(name))+10];
                for j=1:length(lstN)
                    beta=beta-model.Coefficients.Estimate(lstN(j))*modelvals.(model.Coefficients.Row{lstN(j)});
                    UE=UE+m*model.Coefficients.SE(lstN(j));
                    LE=LE-m*model.Coefficients.SE(lstN(j));
                end
                
                lst=find(ismember(model.Coefficients.Row,matlab.lang.makeValidName(lstNonCat{i})));
                lst2=find(abs(modelvals.(model.Coefficients.Row{lst}))>eps(100));
                
                 if(~seperateplot)
                      subplot(mm,m2,1+i);
                     hold on;
                 else

                   figure;
                    set(groot,'defaultTextInterpreter','none')
                    hold on;
                 end

               
               
                if(~isempty(EB))
                    e=nirs.util.scattererrorbar(modelvals.(name)(lst2),modelvals.beta(lst2),EB(lst2));
                    set(e,'linewidth',2,'color','k');
                end     
                     
                s=scatter(modelvals.(name)(lst2),modelvals.beta(lst2));
               
                slope=[model.Coefficients.Estimate(lst) ...
                    model.Coefficients.Estimate(lst)-model.Coefficients.SE(lst)...
                    model.Coefficients.Estimate(lst)+model.Coefficients.SE(lst)];
                XX=[];
                XX(:,1)=LE+x*slope(1);
                XX(:,2)=UE+x*slope(1);
                XX(:,3)=UE+x*slope(2);
                XX(:,4)=LE+x*slope(2);
                XX(:,5)=UE+x*slope(3);
                XX(:,6)=LE+x*slope(3);
                
                
               
                l3=plot(x,x*slope(1));
                L4(1)=plot(x,min(XX,[],2));
                L4(2)=plot(x,max(XX,[],2));
                
                set(s,'SizeData',150,'MarkerFaceColor','k','MarkerEdgeColor','k')
                set(l3,'LineWidth',4,'Color','r','LineStyle','-');
                set(L4,'LineWidth',3,'Color','k','LineStyle','--');
                
                
                a=min(modelvals.(name)(lst2));
                b=max(modelvals.(name)(lst2));
                r=(b-a)/10;
                
                set(gca,'xlim',[a-r b+r]);
                xlabel(lstNonCat{i});
                ylabel(model.ResponseName)
                  set(gca,'fontsize',14)
                title(lstNonCat{i});
            end
            
            
            
            
        end
        
        % critical value
        function out = getCritT( obj, s )
            %% getCritT - returns critical t-stat
            %
            % Args:
            %     s - string specifying statistical significance
            %         (e.g. 'p < 0.05' or 'q < 0.1')
            
            s = strtrim( strsplit( s, '<' ) );
            
            if s{1} == 'p'
                pcrit = str2num(s{2});
                out = - tinv( abs(pcrit)/2, obj.dfe );
                
            elseif s{1} == 'q'
                t = abs(obj.tstat(:))';
                q = obj.q(:);
                
                [~,idx] = unique(q);
                
                out = interp1(q(idx), t(idx), str2num(s{2}));
            end
        end
        
        function out = table( obj )

            if(length(obj)>1)
                out=[];
                for i=1:length(obj)
                    s=obj(i).table;
                    d=nirs.createDemographicsTable(obj(i));
                    d.scanIdx={['scan-' num2str(i)]};
                    d=repmat(d,height(s),1);
                    s=[d s];
                    out=nirs.util.safe_table_vcat({out,s});
                end
                return;
            end


            %% table - returns a table of the regression stats
            beta = obj.beta;
            tstat = obj.tstat;
            p = obj.p;
            q = obj.q;
            se = sqrt(diag(obj.covb));
            dfe = obj.dfe .* ones(size(p));
            
            [minDiscoverableChange,power] = nirs.math.MDC(obj,.8,.05);
            RelativePower=min(minDiscoverableChange)./minDiscoverableChange;
            variab=obj.variables;
            
            %             if(isa(obj.probe,'nirs.core.ProbeROI'))
            %                 [~,i]=ismember(variab(:,1:2),obj.probe.link);
            %                 for j=1:height(obj.probe.link)
            %                     if(iscellstr(obj.probe.link.type{j}) | ischar(obj.probe.link.type{j}))
            %                     N{j}=[obj.probe.link.ROI{j} ':' obj.probe.link.type{j}];
            %                     else
            %                          N{j}=[obj.probe.link.ROI{j} ':' num2str(obj.probe.link.type(j))];
            %                     end
            %                 end
            %
            %                 variab=[table({N{i}}','VariableNames',{'Region'}) variab];
            %                 variab.source=[];
            %                 variab.detector=[];
            %             end
            out = [variab table(beta, se, tstat, dfe, p, q,minDiscoverableChange,RelativePower)];
        end
        
        function out = sorted( obj, colsToSortBy )
            %% sorted - returns sorted stats by columns in variables
            out = obj;
            if nargin < 2
                if( isa(obj(1).probe,'nirs.core.ProbeROI'))
                    colsToSortBy = {'ROI', 'type', 'cond'};
                else
                    colsToSortBy = {'source', 'detector', 'type', 'cond'};
                end
            end
            
            if(length(obj)>1)
                for idx=1:length(obj)
                    out(idx)=sorted(obj(idx),colsToSortBy);
                end
                return
            end
            [out.variables, idx] = nirs.util.sortrows(out.variables, colsToSortBy);
            out.probe.link = nirs.util.sortrows(out.probe.link,{colsToSortBy{ismember(colsToSortBy,out.probe.link.Properties.VariableNames)}}); %out.probe.link(idx,:);
            out.beta = obj.beta(idx);
            out.covb = obj.covb(idx, idx);
            
            if(~isempty(out.pvalue_fixed))
                out.pvalue_fixed=out.pvalue_fixed(idx);
            end
            
        end
        
        [stats,haserror] = ttest( obj, c, b, names );
        stats = ftest( obj, m );
        stats = jointTest( obj );
        
        f = draw( obj, vtype, vrange, thresh,fhandle,ConditionsShown,TypesShown );
        
        printAll( obj, vtype, vrange, thresh, folder, ext );
        
        f = gui(obj);
        
    end
    
    
    methods (Access = protected)
        newNames = transformNames( obj, T );
    end
    
end