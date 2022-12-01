classdef QT < nirs.modules.AbstractModule
%% QT - Performs quality control of the data.
% 
%
% List of parameters:
% 
% Parameter keyword Description
%  freqCut:   1x2 array [fmin fmax] representing the bandpass of the cardiac pulsation (default [0.5 2.5])
%  window :   length in seconds of the window to partition the signal with (defaut: 5)
%  overlap:     fraction overlap (0..0.99) between adjacent windows (default: 0, no overlap)
%  q_threshold:   The required quality value (normalized; 0-1) of good-quality windows in every channel (default: 0.75)
%  conditionsMask:   A binary mask (or the keyword 'all') to indicate the conditions for computing the periods of interest (default: 'all')
%  lambdaMask:    A binary array mapping the selected two wavelength to compute the SCI (default: [1 1], the first two WLs)
%  dodFlag:     A flag indicating to work from DOD data (default: 0)
%  guiFlag:     A flag indicating whether to start or not the GUI.


    properties

        fCut = [0.5 2.5];
        windowSec = 5;
        windowOverlap = 0;
        qThreshold = 0.75;
        lambda_mask = 'all';
        dodFlag = 0;
        guiFlag = 0;
        condMask = 'resting'
        sciThreshold = 0.8;
        pspThreshold = 0.1;
    end
    
    methods

        function obj = QT( prevJob )
           obj.name = 'Quality Testing (QT)';
           if nargin > 0
               obj.prevJob = prevJob;
           end
           obj.citation='';
           
        end
        
        function S = runThis( obj, data )
            for i = 1:numel(data)
                if(isa(data(i),'nirs.core.Data')) 
                    S(i) = nirs.core.QTNirs();
                    %- sort probe.link
                    [~,qt_sort_idx] = sortrows(data(i).probe.link, {'source', 'detector','type'});
                    data(i).probe.link = data(i).probe.link(qt_sort_idx,:);
                    %- sort data
                    data(i).data = data(i).data(:,qt_sort_idx);
                    S(i).probe = data(i).probe;
              
                    qMats = qtnirs(data(i),...
                    'freqCut',obj.fCut,...
                    'window',obj.windowSec,...
                    'overlap',obj.windowOverlap,....
                    'sciThreshold',obj.sciThreshold,...
                    'pspThreshold',obj.pspThreshold,...
                    'qualityThreshold',obj.qThreshold,...
                    'conditionsMask',obj.condMask,...
                    'dodFlag',obj.dodFlag,...
                    'guiFlag',obj.guiFlag);   
                
                    [~,qt_sort_idx] = sortrows(qMats.good_combo_link,[1,2]); %sorting rows by source(1) and detectors (2) columns
                    qMats.good_combo_link = qMats.good_combo_link(qt_sort_idx,:);
                    qMats.sci_array = qMats.sci_array(qt_sort_idx,:);
                    qMats.power_array = qMats.power_array(qt_sort_idx,:);
                    qMats.combo_array = qMats.combo_array(qt_sort_idx,:);
                    qMats.bad_links =  find(mean(qMats.combo_array,2)<qMats.thresholds.quality); 
                    %qMats.bad_links = data(i).probe.link(idx_badlinks,[1 2]);
                    qMats = rmfield(qMats,{'woi','combo_array_expanded','cardiac_data'});                    
                    S(i).qMats = qMats;
                    fprintf('Scan %i of %i processed.\n',i,numel(data));
                else
                    error('unsupported type');
                end
                
            end
            %obj.reportGroup.Properties.VariableNames = {'file_idx','Bad_Links','Bad_Windows'};
        end     
    end  
end

