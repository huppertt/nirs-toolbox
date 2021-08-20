classdef Data < handle
    %% DATA - Holds Real-time NIRS data.
    % 
    % Properties: 
    %     description  - description of data (e.g. filename)
    %     data         - ntime x nchan array of NIRS data
    %     probe        - a Probe object containing measurement geometry
    %     time         - a vector containing the time
    %     Fm           - modulation frequency in MHz (0 for CW NIRS)
    %     Fs           - sampling frequency
    %     stimulus     - a Dictionary object containing stimulus info
    %                    (e.g. stimulus('tapping') returns a StimulusEvents Object)
    %     demographics - a Dictionary object containing demographics info
    %                    (e.g. demographics('age') returns 28)
    % 
    %     
    %  Methods:
    %     draw - displays the time series data and stimulus timings
    
    properties
        description
        data               % channel time series in columns
        probe;             % object describing geometry
        time               % vector of time points
        Fm = 0             % modulation frequency in MHz: 0 for CW; 110 for ISS
        auxillary = Dictionary();  % to hold generic time series information
        stimulus        = Dictionary();	% struct containing stim vectors (vectors, names, types)
        demographics    = Dictionary();	% table containing demographics (names, values)
        updatefunction = {};
    end
    
    properties( Dependent = true )
        Fs = 0;             % sampling frequency in Hz
    end
    properties (Hidden=true)
        rawdata=[];
        rawprobe=[];
        drawhandles_data;
        drawhandles_stimulus;
    end
    
    methods
        
        function reset(obj)
            obj.rawdata=[];
            obj.probe=obj.rawprobe;
            obj.rawprobe=[];
            obj.stimulus=Dictionary();
            obj.auxillary=Dictionary();
            obj.drawhandles_data=[];
            obj.drawhandles_stimulus=[];
        end
        function adddata(obj,d,t)
            
           
            
            if(isempty(obj.time) && ~isempty(obj.probe))
                obj.rawprobe=obj.probe;
                obj.drawhandles_data = zeros(size(d,2),1);
            end
           
            if(~isempty(obj.updatefunction))
                [d2,t,obj.probe,obj.stimulus]=obj.updatefunction.update(d,t,obj.rawprobe,obj.stimulus);
            else
                d2=d;
            end
            obj.rawdata=[obj.rawdata; d];
            obj.data=[obj.data; d2];
            obj.time=[obj.time; t];
            do_once=false;
            for i=1:length(obj.drawhandles_data)
                if(ishandle(obj.drawhandles_data(i)) && obj.drawhandles_data(i)~=0)
                    set(obj.drawhandles_data(i),...
                        'Ydata',[get(obj.drawhandles_data(i),'Ydata') d2(:,i)],...
                        'Xdata',[get(obj.drawhandles_data(i),'Xdata') t]);   
                    if(~do_once)
                        set(get(obj.drawhandles_data(i),'parent'),'Xlim',[obj.time(1) t]);
                        yy=get(get(obj.drawhandles_data(i),'parent'),'Ylim');
                        ymin=min([yy(:); d2(:)]);
                        ymax=max([yy(:); d2(:)]);
                        set(get(obj.drawhandles_data(i),'parent'),'Ylim',[ymin ymax]);
                        do_once=true;
                    end
                end
           end
        end
        
        function data = raw(obj)
            data=convert(obj);
            data.data=obj.rawdata;
            data.probe=obj.rawprobe;
        end
        
        function addeventStart(obj,name,onset,amp)
            if(nargin<4)
                amp=1;
            end
            if(obj.stimulus.iskey(name))
                stim=obj.stimulus(name);
            else
                stim=nirs.design.StimulusEvents;
                stim.name=name;
            end
            
            if(length(stim.dur)>0 && stim.dur(end)==inf)
                warning('stimulus must be Ended before starting a new onset')
                return;
            end
            
            stim.onset(end+1,1)=onset;
            stim.amp(end+1,1)=1;
            stim.dur(end+1,1)=inf;
            obj.stimulus(name)=stim;
            drawstim(obj);
        end
        
        function addeventEnd(obj,name,offset)
            if(obj.stimulus.iskey(name))
                stim=obj.stimulus(name);
            else
                warning('stimulus must be Started before Ending')
                return;
            end
            if(length(stim.dur)>0 && stim.dur(end)==inf)
                if(offset<stim.onset(end))
                    warning('Endtime must be greater than Starttime');
                    return
                end
                stim.dur(end)=offset-stim.onset(end);
            else
                warning('stimulus must be Started before Ending');
                return
            end
            obj.stimulus(name)=stim;
            drawstim(obj);
        end

        
        function addevent(obj,name,onset,dur,amp)
            if(nargin<4)
                dur=1/obj.Fs;
            end
            if(nargin<5)
                amp=1;
            end
            if(obj.stimulus.iskey(name))
                stim=obj.stimulus(name);
            else
                stim=nirs.design.StimulusEvents;
                stim.name=name;
            end
            
            stim.onset(end+1,1)=onset;
            stim.amp(end+1,1)=1;
            stim.dur(end+1,1)=dur;
            obj.stimulus(name)=stim;
            drawstim(obj);
        end
        
        
        function data = convert(obj)
            data=nirs.core.Data;
            data.probe=obj.probe;
            data.stimulus=obj.stimulus;
            data.demographics=obj.demographics;
            data.time=obj.time;
            data.data=obj.data;
            data.description=obj.description;
            data.auxillary=obj.auxillary;
        end
        
        function drawstim(obj)
            for i=1:length(obj.drawhandles_data)
                if(ishandle(obj.drawhandles_data(i)))
                    axis_handle=get(obj.drawhandles_data(i),'parent');
                    obj.draw([],false,axis_handle);
                    return;
                end
            end
        end
        
        %----------------------
        % Same as nirs.core.Data below here
        
        function obj = Data(probe, Fm, stimulus, demographics, description )
            %% Data - Creates a Data object.
            % 
            % Args:
            %     data         - (optional) ntime x nchan array of NIRS data
            %     time         - (optional) a vector containing the time
            %     probe        - (optional) a Probe object containing measurement geometry
            %     Fm           - (optional) modulation frequency (0 for CW NIRS)
            %     stimulus     - (optional) a Dictionary object containing stimulus info
            %     demographics - (optional) a Dictionary object containing demographics info
            %     description  - (optional) description of data (e.g. filename)

            if nargin > 1, obj.probe        = probe;        end
            if nargin > 2, obj.Fm           = Fm;           end
            if nargin > 3, obj.stimulus     = stimulus;     end
            if nargin > 4, obj.demographics = demographics; end
            if nargin > 5, obj.description  = description;  end
        end
        
        function set.stimulus( obj, stim )
           assert( isa(stim,'Dictionary') )
           obj.stimulus = stim;
        end
        
        function set.demographics( obj, demo )
           assert( isa(demo,'Dictionary') )
           obj.demographics = demo;
        end
        
        function set.auxillary( obj, auxillary )
            assert( isa(auxillary,'Dictionary') )
            obj.auxillary= auxillary;
        end
        
        
        function set.time( obj, time )
           assert( isvector(time) | isempty(time))
           obj.time = time(:);
           
        end
        
        function out = get.Fs( obj )
            if length(obj.time) > 1
                out = 1 / mean(diff( obj.time ));
            else
                out = NaN;
            end
        end
        
        
        function showStimMatrix(obj)
            if(isempty(which('spm_DesRep')))
                warning('this function  requires SPM to be installed in this computer');
                return
            end
            [X,conds] = obj.getStimMatrix;
            
            C = ones(size(obj.time));
            
            for id=1:size(C,2)
                cnames{id}=sprintf('trend_%d',i);
            end
            
            SPM.xX.X      = [X C];
            nscan=length(obj.time);
            
            % Create the SPM structure
            
            SPM.xX.iH     = [];
            SPM.xX.iC     = [1:size(X,2)];
            SPM.xX.iB     = [1:size(C,2)] + size(X,2);
            SPM.xX.iG     = [];
            SPM.xX.name   = {conds{:} cnames{:}};
            
            spm_DesRep('DesMtx',SPM.xX);
            
        end
        
        function [X,conds] = getStimMatrix( obj )
            %% returns a [time x condition] matrix of task timings and a cell array of condition names
            if numel(obj)>1
                X = cell(size(obj));
                conds = cell(size(obj));
                for i = 1:numel(obj)
                    [X{i},conds{i}] = obj(i).getStimMatrix;
                end
                return
            end
            t = obj.time;
            conds = obj.stimulus.keys;
           % X = zeros(length(t),length(conds));  %This is less efficient
           % but allows X to keep the class of the stim vector which allows
           % comidations for ordinal and catagorcial vectors
            for i = 1:length(conds)
                stim = obj.stimulus(conds{i});
                X(:,i) = stim.getStimVector(t);
            end
        end
        
        function out = sorted( obj, colsToSortBy )
            %% returns sorted channels of data by column in probe.link
            out = obj;
            if nargin < 2
                colsToSortBy = {'source', 'detector', 'type'};
            end
            if(length(obj)>1)
                for idx=1:length(obj)
                    out(idx)=sorted(obj(idx),colsToSortBy );
                end
                return
            end
            [out.probe.link, idx] = sortrows(out.probe.link, colsToSortBy);
            out.data = out.data(:,idx);
        end
        
        function varargout=draw( obj, lstChannels,adderr,axis_handle )
            %% draw - Plots the probe geometry.
            % 
            % Args:
            %     lstChannels - list of channels to show
            
            if(nargin<4 || isempty(axis_handle))
                axis_handle=gca;
            end
            
            % get data
            if (nargin == 1)
                lstChannels = 1:size(obj(1).data,2);
            end
            if(nargin<3 || isempty(adderr))
                adderr=false;
            end
          
            if(~isempty(obj.drawhandles_stimulus))
                delete(obj.drawhandles_stimulus);
                obj.drawhandles_stimulus=[];
            end
            
            
            if(isempty(lstChannels))
                % draw a plot, but then turn it off to only show the stim
                % marks
                lstChannels=1;
                showplot=false;
            else
                showplot=true;
                obj.drawhandles_data=[];
            end
            if(isempty(obj(1).data))
                lstChannels=[];
            end
            
            if(length(obj)>1)
                figure;
                a=round(sqrt(length(obj)));
                b=ceil(sqrt(length(obj)));
                for i=1:length(obj)
                    subplot(a,b,i);
                    obj(i).draw(lstChannels,adderr);
                    legend off;
                end
                return
            end
            
            t = obj.time;
            d = obj.data(:,lstChannels);
            
            % get stim vecs
            s = []; k = obj.stimulus.keys;
            for i = 1:length( obj.stimulus.keys )
                try
                    s = [s obj.stimulus.values{i}.getStimVector( t )];
                end
            end
            

            % plots
            nextplot=get(axis_handle,'NextPlot');
            hold(axis_handle,'on');
            
            % data min/max/size
            dmax = max( real(d(:)) );
            dmin = min( real(d(:)) );
            dsize = (dmax-dmin);
            
            % plot stim blocks if available
            if ~isempty(s) 
                if(~isa(s,'ordinal') & ~isa(s,'categorical'))
                    for j=1:size(s,2)
                        s(:,j)=s(:,j)./(ones(size(s(:,j)))*max(s(:,j)));
                    end
                    
                    
                    % min/max of axes
                    %pmin = dmin - 0.2*dsize;
                    %pmax = dmin - 0.05*dsize;
                    
                    % adjust amplitude so stims are visible
                    %s = (pmax-pmin)*s + pmin ;
                    s = s*2E7-1E6;
                    % plot
                    obj.drawhandles_stimulus=plot(axis_handle, t, s, 'LineWidth', 3 );
                    
                    % legend
                    l = legend(axis_handle,k{:});
                    set(l,'Interpreter', 'none');
                %    dmax = max( dmax,max(s(:)));
                %    dmin = min( dmin,min(s(:)));
                %    dsize = (dmax-dmin);
                else
                    m=get(axis_handle,'Position');
                    aa=axes('parent',get(axis_handle,'parent'),'Position',[m(1) m(2) m(3) m(4)/8]);
                    obj.drawhandles_stimulus=plot(aa, t, s, 'LineWidth', 3 );
                    set(aa,'YTickLabel',categories(s));
                    set(axis_handle,'Position',[m(1) m(2)+m(4)/8 m(3) m(4)*7/8]);
                    set(axis_handle,'Xtick',[]);
                end
                
                
                
            end
                % min/max of axes
               	pmin = dmin - 0.1*dsize;
                pmax = dmax + 0.1*dsize;
            
            
            % plot data
            if(~isreal(d) & adderr)
                h=errorbar(axis_handle, t, real(d),imag(d) );
            else
                h=plot(axis_handle, t, real(d) );
            end
            obj.drawhandles_data(lstChannels)=h(end:-1:1);
            xlabel(axis_handle, 'seconds' );
            for i = 1:length(obj.stimulus.keys)
                legend(axis_handle,obj.stimulus.keys)
            end
            
                
            if(isempty(t) || min(t)==max(t))
            % axes limits
             xlim(axis_handle, [0 1] );
            else
                xlim(axis_handle, [min(t) max(t)] );
            end
            if(isempty(pmin))
                ylim(axis_handle,[-1 1]);
            elseif (pmin == pmax)
                
                ylim(axis_handle,pmin + [-1 1]);
            else
                ylim(axis_handle, [pmin pmax] );
            end
            
            set(axis_handle,'NextPlot',nextplot);
            
            if(~showplot)
                set(h,'visible','off');
            end
            
            if(nargout>0)
                varargout{1}=h;
            end
            
        end
        
          function varargout=drawwaterfall( obj, lstChannels,adderr,axis_handle )
            %% draw - Plots the probe geometry.
            % 
            % Args:
            %     lstChannels - list of channels to show
            
            if(nargin<4 || isempty(axis_handle))
                axis_handle=gca;
            end
            
            % get data
            if (nargin == 1)
                lstChannels = 1:size(obj(1).data,2);
                
            end
            if(nargin<3 || isempty(adderr))
                adderr=false;
            end
            
            if(isempty(lstChannels))
                % draw a plot, but then turn it off to only show the stim
                % marks
                lstChannels=1;
                showplot=false;
            else
                showplot=true;
            end
            if(isempty(obj(1).data))
                lstChannels=[];
            end
            
            if(length(obj)>1)
                figure;
                a=round(sqrt(length(obj)));
                b=ceil(sqrt(length(obj)));
                for i=1:length(obj)
                    subplot(a,b,i);
                    obj(i).draw(lstChannels,adderr);
                    legend off;
                end
                return
            end
            
            t = obj.time;
            d = obj.data(:,lstChannels);
            
            % get stim vecs
            s = []; k = obj.stimulus.keys;
            for i = 1:length( obj.stimulus.keys )
                s = [s obj.stimulus.values{i}.getStimVector( t )];
            end
            

            % plots
            hold(axis_handle,'on');
            
            % data min/max/size
            dmax = max( real(d(:)) );
            dmin = min( real(d(:)) );
            dsize = (dmax-dmin);
            set(axis_handle,'view',[  -44.0191   74.3316])
            % plot stim blocks if available
            if ~isempty(s) 
                s=s./(ones(size(s))*max(s(:)));
                % min/max of axes
                pmin = dmin;
                pmax = dmax;

                % adjust amplitude so stims are visible
                s = .2*(pmax-pmin)*s + pmin ;
                
                % plot
                waterfall(1, t, s, 'parent',axis_handle );
                
                % legend
                l = legend(axis_handle,k{:});
                set(l,'Interpreter', 'none');
           
            dsize = (dmax-dmin);
                
                
            end
                % min/max of axes
               	pmin = dmin - 0.1*dsize;
                pmax = dmax + 0.1*dsize;
            
            
            % plot data
            cc=colorcube(size(d,2));
            for ii=1:size(d,2)
                h(ii)=plot3(ones(size(t))*(ii+1),t,d(:,ii));
            end
            ylabel(axis_handle, 'seconds' );
            for i = 1:length(obj.stimulus.keys)
                legend(axis_handle,obj.stimulus.keys)
            end
            
                
         
             
            if(isempty(t) || min(t)==max(t))
            % axes limits
             ylim(axis_handle, [0 1] );
            else
                ylim(axis_handle, [min(t) max(t)] );
            end
            if(isempty(pmin))
                zlim(axis_handle,[-1 1]);
            elseif (pmin == pmax)
                 zlim(axis_handle,pmin + [-1 1]);
            else
                zlim(axis_handle, [dmin dmax] );
            end
            set(axis_handle,'xlim',[0 (size(d,2)+2)]);
            set(axis_handle,'YDir','reverse');
            
            names{1}='Stim';
            for i=1:length(lstChannels)
               type=obj.probe.link.type(lstChannels(i));
               if(~iscellstr(type))
                   type=[num2str(type) 'nm'];
               end
                names{i+1}=['Src' num2str(obj.probe.link.source(lstChannels(i))) '-Det' ...
                    num2str(obj.probe.link.detector(lstChannels(i))) ' (' type ')'];
            end
            set(axis_handle,'XTick',[1:size(d,2)+1]);
            set(axis_handle,'XTickLabel',names);
            set(axis_handle,'XTickLabelRotation',300);
            
            if(~showplot)
                set(h,'visible','off');
            end
            
            if(nargout>0)
                varargout{1}=h;
            end
            
          end
        
          f = gui( obj);
        
    end
end

