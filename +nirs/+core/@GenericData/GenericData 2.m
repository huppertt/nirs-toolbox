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
        time               % vector of time points
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

