classdef QTNirs 
%% QTNirs - Class
% 
%
% List of parameters:
% List of methods:


    properties
        qMats = struct;
        probe
    end
    
    methods
        % chPruning(FileStats)
        % Take the 1st-level statistical results and performs a channel
        % pruning process over every scan results. The channel pruning is
        % done by setting the 'beta' values to NaN values to low-quality
        % channels (See qtMats.MeaListAct)
        % Ted posted on the fNIRS Analysis Facebook Group (Nov 16 2020):
        % I just edited (and will update shortly) a bunch of the code to allow channels to be removed from the data. If you replace the channel in the raw.data field with NaN's, then the code should handle it properly and the resulting Stats variable should have a NaN for the beta and tstat for that channel. This works for the OLS, AR-IRLS, and NIRS-SPM models. For the OLS and AR-IRLS, you can now also mask out time points (e.g. mask in both channel or time) using NaN. The time mask does not need to be the same for each channel. Currently you cannot do this with the NIRS-SPM module however because it requires rewriting some of their code to handle this.
        % I also modified the MixedEffects modules (I have not gotten around to the other group level models yet though) that allows NaN in the stats variables. So if subject 1 is missing channel A and subject 2 is missing channel B, etc, then the group model will still average all the valid data together (just the degrees-of-freedom per channel is different).
        % I need to test this all a bit because it was a rather major change to allow this, but I will upload later today.
        % This is also consistent with a commit made in Nov 19 2020 on
        % Github "Added ability to handle NaN in the GLM 1st and 2nd level
        % models" https://github.com/huppertt/nirs-toolbox/commit/4bef48688f116bb537fbb298cf822b0816ae6569
        
        function FileStats = chPruning(obj,FileStats)
            if(nargin<2 || ~isa(FileStats(1),'nirs.core.ChannelStats'))
                error('A ChannelStats array is needed.');
            end
            if (length(obj)~=length(FileStats))
                error('The number of scan-quality and ChannelStats objects do not match.');
            end
            for i = 1:length(FileStats)
                FileStats(i) = FileStats(i).sorted({'source','detector','type'});
                channelsHQ = obj(i).channels2go();                
                idxBadLinks = find(channelsHQ.include==0);
                for j=1:length(idxBadLinks)
                  src = channelsHQ.source(idxBadLinks(j));
                  det = channelsHQ.detector(idxBadLinks(j));
                  vars = FileStats(i).variables;
                  idxVars = vars.source == src & vars.detector ==det;
                  FileStats(i).beta(idxVars) = NaN;
%                  FileStats(i).tstat(idxVars) = NaN;                  
                end
            end         
        end
        
        
        % scans2go(min_link_perc)
        % returns a mask with those scans with at least min_link_perc% of
        % channels with good quality. A channel is considered good if
        % its quality is greater than the 'quality threshold'. 
        % For instance the quality of a channel 'i' is evaluated as:
        % qMats.good_combo_link(i,3) >= qMats.thresholds.quality        
        function scansHQ = scans2go(obj,min_link_perc)
            if(nargin<2 || isempty(min_link_perc)), min_link_perc = 0.9; end
            if length(obj)> 1
                scansHQ = false(length(obj),1);
                for i =1: length(obj)
                    thrQty = obj(i).qMats.thresholds.quality;
                    scansHQ(i) = sum(obj(i).qMats.good_combo_link(:,3)>=thrQty)/size(obj(i).qMats.good_combo_link,1) >= min_link_perc;
                end
            end
        end
        
        % channels2go()
        % returns a table with surviving channels after channel prunning 
        % process based on quality assesment.
        % obj.qMats.MeasListAct contains a mask with good-quality channels
        function channelsHQ = channels2go(obj)           
            for i =1: length(obj)    
                %channelsHQ{i} = table();
                channelsHQ{i} = array2table([obj(i).qMats.MeasList(:,[1 2 4]),...
                    obj(i).qMats.MeasListAct]);
                channelsHQ{i}.Properties.VariableNames={'source','detector','type','include'};
                [channelsHQ{i}, idx] = sortrows(channelsHQ{i},{'source','detector','type'});               
            end
            if length(obj)==1
                channelsHQ = channelsHQ{1};
            end
        end
        
        
        function f=drawBarBadLinks(obj)            
            report_mat=zeros(1,size(obj(1).probe.link,1)/2);
            master_threshold = obj(1).qMats.thresholds.quality;
            for i =1:length(obj)
               report_mat = report_mat + (obj(i).qMats.good_combo_link(:,3)>=master_threshold)';
            end
            %report_mat = report_mat./max(report_mat);
            f=figure('Name','Channel-level Report','NumberTitle','off');
            bar(report_mat);
            title('High-quality Channels');
            xticks(1:size(obj(1).probe.link,1)/2);
            xticklabels(string([num2str(obj(1).qMats.good_combo_link(:,1)),...
                repmat('-',size(obj(1).qMats.good_combo_link,1),1),...
                num2str(obj(1).qMats.good_combo_link(:,2),2)]));
            xtickangle(35);
            xlabel('Channel #');
            ylabel('#Scans');
            legend(['HQ above ',num2str(master_threshold*100),'% of time']);
        end
        
        function f=draw(obj,type)
            % types can be: SNI, SNR, Mean, Median, Max, Min, Motion, KPSS,
            % ADtest
            
            types={'sq','sci','psp'};
            
            if(nargin<2)
                type='sq';
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
                case 'sci'
                    val = mean(obj.qMats.sci_array,2);
                    vrange=[min(val) max(val)];
                    idx_below_thresh = val<obj.qMats.thresholds.sci;
                case 'psp'
                    val = mean(obj.qMats.power_array,2);
                    vrange=[min(val) max(val)];
                    idx_below_thresh = val<obj.qMats.thresholds.peakpower;
                case 'sq'
                    val = obj.qMats.good_combo_link(:,3);
                    vrange=[0 1];
                    idx_below_thresh = val<obj.qMats.thresholds.quality;
            end
            
            %vrange=[min(val) max(val)];
           
            [~,cmap] = evalc(' cbrewer(''seq'',''Greens'',128) ');
            z = linspace(vrange(1), vrange(2), size(cmap,1))';
            %idx_below_thresh = z<obj.qMats.thresholds.quality; 
            %cmap(idx_below_thresh,:) = repmat(cmap(43,:),sum(idx_below_thresh),1);
             
            lstyles = repmat({'LineStyle', '-', 'LineWidth', 8}, [length(val) 1]);
            lstyles(idx_below_thresh,:) = repmat({'LineStyle', '--', 'LineWidth', 4},[sum(idx_below_thresh),1]);
            idx = bsxfun(@minus, val', z);
            [~, idx] = min(abs(idx), [], 1);
            %idx_badchannels = find(z<reportTable.thresholds.quality);
            %cmap(idx_badchannels,:) = repmat(cmap(1,:),length(idx_badchannels),1);
            colors = cmap(idx, :);
            
            u=unique(obj.probe.link.type);
            f=figure;
            
            %lst=ismember(obj.probe.link.type,u(1));
            h=axes;%subplot(length(u),1,i);
            %obj.probe.draw(colors(lst,:),lstyles,h);
            obj.probe.draw(colors,lstyles,h);
            c = colorbar; colormap(cmap); caxis(vrange);
            if strcmp(type,'sq') || strcmp(type,'sci')
                c.TickLabels = split(num2str(c.Ticks*100));
            else
                c.TickLabels = split(num2str(c.Ticks));
            end
            c.Label.String='(%)';
            
            title(['Channel quality:', type]);
        end
        
       function f=drawGroup(obj,type)
            % types can be: 'sqmask','sqcount','sq','sci','psp','bar'
            % 
            
            % -- Modified by SAMH
            types={'sqmask','sqcount','sq','sci','psp','bar',};
            
            if(nargin<2)
                type='sq';
            end
            if(~ismember(lower(type),lower(types)))
                disp('type must be one of:');
                disp(types);
                return
            end
            
            if(length(obj)>1)
                gcl_mean = zeros(size(obj(1).qMats.good_combo_link,1),1);
                sci_mean = zeros(size(obj(1).qMats.good_combo_link,1),1);
                psp_mean = zeros(size(obj(1).qMats.good_combo_link,1),1);
                reportGroup = table;
                
                for i=1:length(obj)
                    gcl_mean = gcl_mean + obj(i).qMats.good_combo_link(:,3);
                    sci_mean = sci_mean + mean(obj(i).qMats.sci_array,2);
                    psp_mean = psp_mean + mean(obj(i).qMats.power_array,2);
                    new_table_row = table(i, {obj(i).qMats.bad_links'}, {obj(i).qMats.bad_windows});
                    reportGroup = [reportGroup; new_table_row];                                      
                end
                gcl_mean = gcl_mean /length(obj);
                sci_mean = sci_mean /length(obj);
                psp_mean = psp_mean /length(obj);
                reportGroup.Properties.VariableNames = {'scan_idx','Bad_Links','Bad_Windows'};                      
            elseif (length(obj)==1)
                obj.draw(type);
                return;
            end
            
            switch(lower(type))
                case 'sci'
                    val = sci_mean;
                    vrange=[min(val) max(val)];
                    idx_below_thresh = val<obj(1).qMats.thresholds.sci;
                    
                case 'psp'
                    val = psp_mean;
                    vrange=[min(val) max(val)];
                    idx_below_thresh = val<obj(1).qMats.thresholds.peakpower;
                    
                case 'bar'
                    f=obj.drawBarBadLinks();
                    return;
                case 'sq'
                    val = gcl_mean;
                    vrange=[0 1];
                    idx_below_thresh = val<obj(1).qMats.thresholds.quality;
                
                % FIX THIS VISUALIZATION ------------------------------------
                case 'sqcount'
                    disp('Not implemented, please use other option.');
                    return;
%                     % combo results of the whole group in percentage (0-bad
%                     % windows, 1-good windows)
%                     % Counting at window level
%                     sq_count = zeros(size(obj(1).qMats.combo_array));
%                     for i = 1:length(obj)
%                         sq_count = sq_count + obj(i).qMats.combo_array;
%                     end
%                     val = sum(sq_count,2)/(size(sq_count,2)*length(obj))*100;
%                     %vrange=[min(val) max(val)];
%                     vrange=[0 100];
%                     idx_below_thresh = false(length(val),1);
                %------------------------------------------------------------
                case 'sqmask'                    
                    % Mask of channels above the threshold, for every
                    % subject, a channel quality above the threshold (90%)
                    % adds '1' to the total. The sum is normalized by the number
                    % of scans, the expressed as percentage.
                    % Counting at Channel level
                    gcl_mask = zeros(size(obj(1).qMats.good_combo_link,1),1);
                    for i=1:length(obj)
                        gcl_mask = gcl_mask + (obj(i).qMats.good_combo_link(:,3)>=obj(1).qMats.thresholds.quality);
                    end
                    %gcl_mask = gcl_mask/length(obj)*100;
                    gcl_mask = gcl_mask;
                    val = gcl_mask;
                    %vrange=[min(val) max(val)];
                    vrange=[0 length(obj)];
                    idx_below_thresh = false(length(val),1);
            end
            
                        %vrange=[min(val) max(val)];
            [~,cmap] = evalc(' cbrewer(''div'',''RdYlGn'',128) ');
            z = linspace(vrange(1), vrange(2), size(cmap,1))';
            %idx_below_thresh = z<obj.qMats.thresholds.quality; 
            %cmap(idx_below_thresh,:) = repmat(cmap(43,:),sum(idx_below_thresh),1);
            lstyles = repmat({'LineStyle', '-', 'LineWidth', 8}, [length(val) 1]);
            lstyles(idx_below_thresh,:) = repmat({'LineStyle', '--', 'LineWidth', 4},[sum(idx_below_thresh),1]);
            idx = bsxfun(@minus, val', z);
            [~, idx] = min(abs(idx), [], 1);
            %idx_badchannels = find(z<reportTable.thresholds.quality);
            %cmap(idx_badchannels,:) = repmat(cmap(1,:),length(idx_badchannels),1);
            colors = cmap(idx, :);
            
            u=unique(obj(1).probe.link.type);
            f=figure;
            
            %lst=ismember(obj(1).probe.link.type,u(1));
            h=axes;%subplot(length(u),1,i);
            %obj(1).probe.draw(colors(lst,:),lstyles,h);
            obj(1).probe.draw(colors,lstyles,h);
            c = colorbar; colormap(cmap); caxis(vrange);
           if strcmp(type,'sq')
                c.TickLabels = split(num2str(c.Ticks*100));
                c.Label.String='High-quality scans (%)';
            else
                c.TickLabels = split(num2str(c.Ticks));
                c.Label.String='High-quality scans';
            end
            
            c.Label.FontSize = 12;
            
            title(['Channel quality:', type]);
        end        
    end
    
end

