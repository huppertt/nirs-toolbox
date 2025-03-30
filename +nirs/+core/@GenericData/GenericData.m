classdef GenericData
    %% DATA - Holds generic time series data (eg auxillary data).
    % 
    % Properties: 
    %     description  - description of data (e.g. filename)
    %     data         - ntime x nchan array of NIRS data
    %     link        -  data type info
    %     time         - a vector containing the time
    %     Fs           - sampling frequency
    % 
    %     
    %  Methods:
    %     draw - displays the time series data and stimulus timings
    
    properties
        description
        data               % channel time series in columns
        link;             % object describing geometry
        time;               % vector of time points
        stimulus;
      end
    
    properties( Dependent = true )
        Fs = 0;             % sampling frequency in Hz
    end

    methods
        function obj = GenericData( data, time,link, description )
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

            if nargin > 0, obj.data         = data;         end
            if nargin > 1, obj.time         = time;         end
            if nargin > 2, obj.link        = link;        end
            if nargin > 3, obj.description  = description;  end
            obj.stimulus=Dictionary;

        end
        
        
        function obj = set.time( obj, time )
           assert( isvector(time) )
           obj.time = time(:);
        end
        
        function out = get.Fs( obj )
            if length(obj.time) > 1
                out = 1 / mean(diff( obj.time ));
            else
                out = NaN;
            end
        end
                
        function out = sorted( obj, colsToSortBy )
            %% returns sorted channels of data by column in probe.link
            out = obj;
            if nargin < 2
                colsToSortBy = {'name', 'type'};
            end
            if(length(obj)>1)
                for idx=1:length(obj)
                    out(idx)=sorted(obj(idx),colsToSortBy );
                end
                return
            end
            [out.link, idx] = sortrows(out.link, colsToSortBy);
            out.data = out.data(:,idx);
        end
        
        function varargout=draw( obj, lstChannels,adderr )
            %% draw - Plots the probe geometry.
            % 
            % Args:
            %     lstChannels - list of channels to show
            
            % get data
            if nargin == 1
                lstChannels = 1:size(obj(1).data,2);
                
            end
            if(nargin<3)
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
            
            % plots
            gca; hold on;
            
            % data min/max/size
            dmax = max( real(d(:)) );
            dmin = min( real(d(:)) );
            dsize = (dmax-dmin);
            legendNames={};

            % plot stim blocks if available
            s = []; k = obj.stimulus.keys;
            for i = 1:length( obj.stimulus.keys )
                try
                    s = [s obj.stimulus.values{i}.getStimVector( t )];
                    legendNames{end+1}=obj.stimulus.keys{i};
                end
            end
            if ~isempty(s)
                axis_handle=gca;
                if(~isa(s,'ordinal') & ~isa(s,'categorical'))
                    for j=1:size(s,2)
                        s(:,j)=s(:,j)./(ones(size(s(:,j)))*max(s(:,j)));
                    end


                    % min/max of axes
                    pmin = dmin - 0.2*dsize;
                    pmax = dmin - 0.05*dsize;

                    % adjust amplitude so stims are visible
                    s = (pmax-pmin)*s + pmin ;
                    % plot
                    plot(axis_handle, t, s, 'LineWidth', 3 );

                    % legend
                    l = legend(axis_handle,k{:});
                    set(l,'Interpreter', 'none');
                    dmax = max( dmax,max(s(:)));
                    dmin = min( dmin,min(s(:)));
                    dsize = (dmax-dmin);
                else
                    m=get(axis_handle,'Position');
                    aa=axes('parent',get(axis_handle,'parent'),'Position',[m(1) m(2) m(3) m(4)/8]);
                    plot(aa, t, s, 'LineWidth', 3 );
                    set(aa,'YTickLabel',categories(s));
                    set(axis_handle,'Position',[m(1) m(2)+m(4)/8 m(3) m(4)*7/8]);
                    set(axis_handle,'Xtick',[]);
                end
           end
            
           if(~isempty(obj.link))
               for i=1:length(lstChannels)
                   legendNames{end+1}=obj.link.type{lstChannels(i)};
               end
           end
            % min/max of axes
            pmin = dmin - 0.1*dsize;
            pmax = dmax + 0.1*dsize;

            % plot data
            if(~isreal(d) & adderr)
                h=errorbar( t, real(d),imag(d) );
            else
                h=plot( t, real(d) );
            end
            xlabel( 'seconds' );

            if(~isempty(obj.link))
                if(ismember('units',obj.link.Properties.VariableNames))
                    ylabel(obj.link.units{1});
                end
                if(ismember('type',obj.link.Properties.VariableNames))
                    legend(legendNames);
                end
            end
            
            
            
            % axes limits
            xlim( [min(t) max(t)] );
            if pmin == pmax
                ylim(pmin + [-1 1]);
            else
                ylim( [pmin pmax] );
            end
            
            if(~showplot)
                set(h,'visible','off');
            end
            
            if(nargout>0)
                varargout{1}=h;
            end
            
        end
        
    end
end

