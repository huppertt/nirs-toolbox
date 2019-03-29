classdef ProbeROI < nirs.core.Probe
    %% ProbeROI.  This is a class derived from the core Probe class which adds
    % methods for registration, display, and co-registration based on the 10-20
    % landmarking
    
    properties
        RegionNames
        defaultdrawfcn; % Drawing function to use
    end
    
    methods
        
        function obj = ProbeROI(regionnames)
            
            if(nargin==1)
                obj.RegionNames=regionnames;
            end
            
        end
        
         function obj = set.defaultdrawfcn(obj,str)
            
            if(nargin==1)
                return
            end
            
            allowed={'Patch','Color coded boxes';...
                'Bar', 'Bar Chart'};
            
            if(~isempty(str))
                idx=find(ismember(lower({allowed{:,1}}),lower(str)));
            else
                idx=[];
            end
            
            if(~isempty(idx))
                obj.defaultdrawfcn=allowed{idx,1};
            elseif(isempty(str) || strcmp(str,'?'))
                disp('Here are the options for drawing configurations');
                disp(allowed);
            end
        end
        
        function obj = set.RegionNames(obj,regionnames)
            
            
            for i=1:length(regionnames)
                l=max(strfind(regionnames{i},':'));
                if(~isempty(l))
                    type{i,1}=regionnames{i}(l+1:end);
                    regions{i,1}=regionnames{i}(1:l-1);
                else
                    type{i,1}='unknown';
                    regions{i,1}=regionnames{i};
                end
            end
            
            numsrc=length(unique(regions));
            ur=unique(regions);
            for i=1:length(regions)
                src(i,1)=find(ismember(ur,regions{i}));
            end
            
            %Create a dummy optodes table
            
            warning('off','MATLAB:table:RowsAddedExistingVars');
            for i=1:numsrc
                s=['0000' num2str(i)];
                obj.optodes.Name{i,1}=['Src-' s(end-3:end)];
                obj.optodes.X(i,1)=NaN;
                obj.optodes.Y(i,1)=NaN;
                obj.optodes.Z(i,1)=NaN;
                obj.optodes.Type{i,1}='Source';
                obj.optodes.Units{i,1}='mm';
            end
            obj.optodes.Name{numsrc+1,1}='Det-0001';
            obj.optodes.X(numsrc+1,1)=NaN;
            obj.optodes.Y(numsrc+1,1)=NaN;
            obj.optodes.Z(numsrc+1,1)=NaN;
            obj.optodes.Type{numsrc+1,1}='Detector';
            obj.optodes.Units{numsrc+1,1}='mm';
            
            
            obj.link=table(regions,type,...
                'VariableNames',{'ROI','type'});
            obj.RegionNames=ur;
            
            obj.fixeddistances=ones(height(obj.link),1);
            
        end
        
        function varargout=draw( obj, colors, lineStyles, axis_handle,conditions,vals )
             if(nargin<5)
                conditions={''};
             end
            conditions=unique(conditions);
            n = length(obj.RegionNames)*length(conditions);
            if nargin < 2 || isempty(colors)
                colors = repmat([0.3 0.5 1], [n 1]);
            elseif size(colors,1) == 1
                colors = repmat(colors, [n 1]);
            end  
            if nargin < 3 || isempty(lineStyles)
               lineStyles = repmat({'LineStyle', '-', 'LineWidth', 2}, [n 1]);
           elseif size(lineStyles, 1) < n
               lineStyles = repmat({'LineStyle', '-', 'LineWidth', 2}, [n 1]);
           end
           
           if nargin < 4
               axis_handle = axes();
           end
           
           
            if(nargin<6)
                vals=ones(n,1);
            end
            
            
            if(strcmp(lower(obj.defaultdrawfcn),'patch'))
                varargout= obj.drawPatch(colors, lineStyles, axis_handle,conditions );
            else
                varargout= obj.drawBar(vals, lineStyles, axis_handle,conditions);
            
            end
            if(~iscell(varargout))
                varargout={varargout};
            end
            
            
        end
        
        
       function varargout=drawBar( obj, vals, lineStyles, axis_handle,conditions )
           n = length(obj.RegionNames)*length(conditions);
           if nargin < 3 || isempty(lineStyles)
               lineStyles = repmat({'LineStyle', '-', 'LineWidth', 2}, [n 1]);
           elseif size(lineStyles, 1) < n
               lineStyles = repmat({'LineStyle', '-', 'LineWidth', 2}, [n 1]);
           end
           
           if nargin < 4
               axis_handle = axes();
           end
           vals=reshape(vals,length(obj.RegionNames),length(conditions));
           hold(axis_handle,'on');
           for i=1:size(vals,1)
               h(i)=bar(axis_handle,i,vals(i,:),'grouped');
           end
           set(gca,'XtickLabel',[]);
           legend(obj.RegionNames);
            if(nargout==1)
                varargout{1}=h;
            end
           
       end
        
        
        function varargout=drawPatch( obj, colors, lineStyles, axis_handle,conditions )
            %% draw - Plots the probe geometry.
            %
            % Args:
            %     colors      - (optional) n x 3 array of colors [R, G, B] for each channel
            %     lineStyles  - (optional) 2D cell array with each row containing
            %                   arguments for the 'line' functions (e.g. {'LineWidth',6})
            %     axis_handle - (optional) handle to axis to the plot to
            
            
            n = length(obj.RegionNames)*length(conditions);
            
            if nargin < 2 || isempty(colors)
                colors = repmat([0.3 0.5 1], [n 1]);
            elseif size(colors,1) == 1
                colors = repmat(colors, [n 1]);
            end
            
            if nargin < 3 || isempty(lineStyles)
                lineStyles = repmat({'LineStyle', '-', 'LineWidth', 2}, [n 1]);
            elseif size(lineStyles, 1) < n 
                lineStyles = repmat({'LineStyle', '-', 'LineWidth', 2}, [n 1]);
            end
            
            if nargin < 4
                axis_handle = axes();
            end
            
           
            hold on;
            for i=1:n
                h(i)=fill([i-.2 i-.2 i+.2 i+.2],...
                    [0 1 1 0],colors(i,:));
                 set(h(i),lineStyles{i,:});
                 if(strcmp(get(h(i),'LineStyle'),'--'))
                     set(h(i),'FaceAlpha',.1);
                 end
            end
            set(axis_handle,'Xlim',[0 n+1],'Ylim',[-.1 1.1])
            set(axis_handle,'Xtick',[1:n]);
            set(axis_handle,'Ytick',[]);
            
            nm={}; cnt=1;
            for ii=1:length(conditions)
                for jj=1:length(obj.RegionNames)
                    nm{cnt}=[obj.RegionNames{jj} ':' conditions{ii}];
                    cnt=cnt+1;
                end
            end
            
            set(gca,'XtickLabel',nm);
            
            if(nargout==1)
                varargout{1}=h;
            end
            
        end
        
    end
    
end