classdef QCreport
    %     % QCreport holds information about the quality control of a data object
    %
    %     Properties
    %         description   - description of data (copied from data)
    %         demographics  - copied from data
    %         probe - copied from data
    %
    %     Methods
    %         draw
    %         table
    %         disp
    
    properties
        description;
        probe;
        demographics;
        
    end
    properties(Hidden = true)
        data
        class
        time
        stimulus
        inn
    end
    
    methods
        function obj = QCreport(data)
            if(length(data)>1)
                for i=1:length(data)
                    obj(i)=nirs.core.QCreport(data(i));
                end
                return
            end
            
            obj.description=data.description;
            obj.demographics=data.demographics;
            obj.probe=data.probe;
            
            obj.data=data.data;
            obj.stimulus=data.stimulus;
            obj.class=class(data);
            obj.time=data.time;
            obj.inn = nirs.math.innovations(obj.data,12);
        end
        
        function disp(obj)
            
            if(length(obj)>1)
                disp(['nirs.core.QCreport<' num2str(size(obj,1)) ','...
                    num2str(size(obj,2)) '>']);
                return;
            end
            
            disp(['NIRS Quality Report (' obj.class ')']);
            disp('------------------------------------');
            if(~isempty(obj.description))
                disp(['Description: ' obj.description]);
            end
            disp(['Scan length: ' num2str(obj.time(end)-obj.time(1)) 's '...
                ' [ ' num2str(length(obj.time)) ' samples ]']);
            u=unique(obj.probe.link.type);
            str='';
            for i=1:length(u)
                if(iscellstr(u))
                    str=[str ' ' u{i}];
                else
                    str=[str ' ' num2str(u(i)) 'nm'];
                end
            end
            disp(['Number of Channels: ' num2str(height(obj.probe.link)) ...
                ' [' str ' ]']);
            if(strcmp(obj.class,'nirs.core.Data') & isnumeric(u))
                SNR=median(obj.data)./std(obj.data);
                disp(['SNR: ' num2str(median(SNR)) ' [ ' num2str(min(SNR)) '-'...
                    num2str(max(SNR)) ' ]']);
                for i=1:length(u)
                    if(iscellstr(u))
                        str=u{i};
                    else
                        str=[num2str(u(i)) 'nm'];
                    end
                    lst=find(ismember(obj.probe.link.type,u(i)));
                    SNRa=median(obj.data(:,lst))./std(obj.data(:,lst));
                    disp(['     ' num2str(median(SNRa)) ' [ ' num2str(min(SNRa)) '-'...
                        num2str(max(SNRa)) ' ] @' str]);
                end
                tcric=2*-tinv(0.05,size(obj.data,2));
                disp(['  ' num2str(length(find(zscore(SNR)>tcric))) ' outliers @p<0.05']);
                
            end
            
            %-------
            SNI=mad(obj.data,1,1)./mad(obj.inn,1,1);
            disp(['SNI: ' num2str(median(SNI)) ' [ ' num2str(min(SNI)) '-'...
                num2str(max(SNI)) ' ]']);
            for i=1:length(u)
                if(iscellstr(u))
                    str=u{i};
                else
                    str=[num2str(u(i)) 'nm'];
                end
                lst=find(ismember(obj.probe.link.type,u(i)));
                SNIa=mad(obj.data(:,lst),1,1)./mad(obj.inn(:,lst),1,1);
                disp(['     ' num2str(median(SNIa)) ' [ ' num2str(min(SNIa)) '-'...
                    num2str(max(SNIa)) ' ] @' str]);
            end
            tcric=2*-tinv(0.05,size(obj.data,2));
            disp(['  ' num2str(length(find(zscore(SNI)>tcric))) ' outliers @p<0.05']);
            
            %-------
            tcric=2*-tinv(0.05,length(obj.data));
            mocoevents(1)=length(find(max(zscore(obj.inn)')>tcric))/length(obj.data);
            tcric=2*-tinv(0.01,length(obj.data));
            mocoevents(2)=length(find(max(zscore(obj.inn)')>tcric))/length(obj.data);
            disp(['Motion events: ' num2str(100*mocoevents(1)) '% @ p<0.05']);
            disp(['               ' num2str(100*mocoevents(2)) '% @ p<0.01']);
            
            
        end
        
        function tbl = table(obj)
            u=unique(obj.probe.link.type);
            if(strcmp(obj.class,'nirs.core.Data') & isnumeric(u))
                SNR=(median(obj.data)./std(obj.data))';
            else
                SNR=repmat({'N/A'},height(obj.probe.link),1);
            end
            
            SNI = (mad(obj.data,1,1)./mad(obj.inn,1,1))';
            
            StdErr = std(obj.data)';
            Mean = mean(obj.data)';
            Median = median(obj.data)';
            
            tcric=2*-tinv(0.05,size(obj.data,2));
            Outlier = (zscore(SNI)>tcric);
            
            
            tiny_s = 1e-6 * std(obj.inn);
            t = 4.685;
            s = median(abs(obj.inn)) / 0.6745;
            r=obj.inn./(ones(size(obj.inn,1),1)*(max(s,tiny_s)*t));
            w = (abs(r)<1) .* (1 - r.^2).^2;
            w=sqrt(w);
            MotionFraction=(1-sum(w,1)/length(w))';
            
            for i=1:size(obj.data,2)
                [~,Anderson_Darling(i,1)]=adtest(obj.inn(:,i));
                [~,p]=kpsstest(obj.inn(:,i));
                if(p==0.1)
                    KPSS{i,1}='--';
                else
                    KPSS{i,1}=num2str(p);
                end
            end
            
            tbl=[obj.probe.link table(Mean,Median,StdErr,SNR,SNI,Anderson_Darling,KPSS,MotionFraction,Outlier)];
            
            
        end
        
        function draw(obj,type)
            % types can be: SNI, SNR, Mean, Median, Max, Min, Motion, KPSS,
            % ADtest
            
            types={'SNI','SNR','Mean','Median','Max','Min','Motion','KPSS','ADtest'};
            if(nargin<2)
                type='SNR';
            end
            if(~ismember(lower(type),lower(types)))
                disp('type must be one of:');
                disp(types);
                return
            end
            
            if(length(obj)>1)
                for i=1:length(obj)
                    obj(i).draw(type);
                end
                return
            end
            
            switch(lower(type))
                case 'sni'
                    val=(mad(obj.data,1,1)./mad(obj.inn,1,1))';
                case 'snr'
                    val=(median(obj.data)./std(obj.data))';
                case 'mean'
                    val=(mean(obj.data))';
                case 'median'
                    val=(median(obj.data))';
                case 'max'
                    val=max(obj.data,[],1)';
                case 'min'
                    val=min(obj.data,[],1)';
                case 'motion'
                    tcric=2*-tinv(0.05,size(obj.data,2));
                    SNI=(mad(obj.data,1,1)./mad(obj.inn,1,1))';
                    
                    tiny_s = 1e-6 * std(obj.inn);
                    t = 4.685;
                    s = median(abs(obj.inn)) / 0.6745;
                    r=obj.inn./(ones(size(obj.inn,1),1)*(max(s,tiny_s)*t));
                    w = (abs(r)<1) .* (1 - r.^2).^2;
                    w=sqrt(w);
                    val=(1-sum(w,1)/length(w))';
                case 'kpss'
                    for i=1:size(obj.data,2)
                        [~,val(i,1)]=kpsstest(obj.inn(:,i));
                    end
                case 'adtest'
                    for i=1:size(obj.data,2)
                        [~,val(i,1)]=adtest(obj.inn(:,i));
                    end
                    
            end
            
            vrange=[min(val) max(val)];
            [~,cmap] = evalc('flipud( cbrewer(''div'',''RdBu'',128) )');
            z = linspace(vrange(1), vrange(2), size(cmap,1))';
            idx = bsxfun(@minus, val', z);
            [~, idx] = min(abs(idx), [], 1);
            colors = cmap(idx, :);
            
            u=unique(obj.probe.link.type);
            figure;
            for i=1:length(u);
                lst=ismember(obj.probe.link.type,u(i));
                h=subplot(length(u),1,i);
                obj.probe.draw(colors(lst,:),[],h);
                c = colorbar; colormap(cmap); caxis(vrange);
                if(iscellstr(u))
                    str=u{i};
                else
                    str=[num2str(u(i)) 'nm'];
                end
                title([type ':' str]);
            end
            
        end
        
        
        function plot(obj,type,chan)
            
            types={'motion','innovations','autocorr','dist',...
                'logdist','student','integral',...
                'dist-innov',...
                'logdist-innov','student-innov','integral-innov'};
            
            if(nargin<2 || isempty(type))
                type='dist';
            end
            if(nargin<3)
                chan=1:height(obj(1).probe.link);
            end
            
            if(~ismember(lower(type),lower(types)))
                disp('type must be one of:');
                disp(types);
                return
            end
            
            if(length(obj)>1)
                for i=1:length(obj)
                    obj(i).plot(type,chan);
                end
                return
            end
            
            switch(lower(type))
                case 'motion'
                    tiny_s = 1e-6 * std(obj.inn);
                    t = 4.685;
                    s = median(abs(obj.inn)) / 0.6745;
                    r=obj.inn./(ones(size(obj.inn,1),1)*(max(s,tiny_s)*t));
                    w = (abs(r)<1) .* (1 - r.^2).^2;
                    d=sqrt(w);
                    t=obj.time;
                case 'innovations'
                    d=obj.inn;
                    t=obj.time;
                case 'autocorr'
                    for i=1:size(obj.data,2)
                        [d(:,i),t] = autocorr(obj.data(:,i));
                    end
                case 'dist'
                    dd=obj.data;
                    dd=dd-ones(size(dd,1),1)*mean(dd,1);
                    
                    t=linspace(min(dd(:)), max(dd(:)),30);
                    for i=1:size(dd,2)
                        d(:,i)=hist(dd(:,i),t);
                    end
                case 'logdist'
                    dd=obj.data;
                    dd=dd-ones(size(dd,1),1)*mean(dd,1);
                    t=linspace(min(dd(:)), max(dd(:)),30);
                    for i=1:size(dd,2)
                        d(:,i)=hist(dd(:,i),t);
                    end
                    d=log(d);
                case 'dist-innov'
                    
                    dd=obj.inn;
                    dd=dd-ones(size(dd,1),1)*mean(dd,1);
                    
                    t=linspace(min(dd(:)), max(dd(:)),30);
                    for i=1:size(dd,2)
                        d(:,i)=hist(dd(:,i),t);
                    end
                case 'logdist-innov'
                    dd=obj.inn;
                    dd=dd-ones(size(dd,1),1)*mean(dd,1);
                    
                    t=linspace(min(dd(:)), max(dd(:)),30);
                    for i=1:size(dd,2)
                        d(:,i)=hist(dd(:,i),t);
                    end
                    d=log(d);
                case 'student'
                    d=zscore(obj.data);
                    t=obj.time;
                    
                case 'student-innov'
                    d=zscore(obj.inn);
                    t=obj.time;
                case 'integral'
                    d=obj.data;
                    d=d-ones(size(d,1),1)*mean(d,1);
                    d=cumsum(d);
                    t=obj.time;
                case 'integral-innov'
                    d=obj.inn;
                    d=d-ones(size(d,1),1)*mean(d,1);
                    d=cumsum(d);
                    t=obj.time;
            end
            clear dd;
            figure;
            plot(t,d(:,chan),'color',[.6 .6 .6]);
            hold on;
            if(length(chan)>2)
                
                for i=1:size(d,1);
                    dd(i,:)=sort(d(i,chan));
                end;
                lb=dd(:,round(length(chan)/4));
                ub=dd(:,round(length(chan)*3/4));
                plot(t,lb,'r--');
                plot(t,ub,'r--');
                plot(t,median(d(:,chan),2),'k','linewidth',4);
            end
            title(type);
            
            
        end
        
        
        
        
    end
    
end