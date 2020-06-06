classdef ApplyROI < nirs.modules.AbstractModule
    %% ApplyROI - Performs ROI averaging
    % Usage:
    % job = nirs.modules.ApplyROI();
    % job.listOfROIs(1,:) = array2table({[2 2],[1 2],'left lPFC'});
    % job.listOfROIs(2,:) = array2table({[2 2],[1 2],'left aPFC'});
    % job.listOfROIs(3,:) = array2table({[3 3],[4 5],'right aPFC'});
    % job.listOfROIs(4,:) = array2table({[3 4],[6 6],'right lPFC'});
    % dataROI = job.run( dataChannel );
    %
    % Options:
    %     listOfROIs - [#ROI x 3] table of ROIs (source, detector, name)
    
    
    properties
        listOfROIs = table({},{},{},'VariableNames',{'source','detector','Name'});
        weighted=false;
    end
    
    methods
        function obj = ApplyROI( prevJob )
            if nargin > 0, obj.prevJob = prevJob; end
            obj.name = 'Apply ROI';
        end
        
        % Apply ROI averaging to data and update probe
        function dataROI = runThis( obj , dataChannel )
            

            
            
            
            
            for i=1:length(dataChannel)
                if(~istable(obj.listOfROIs))
                    
                    listOfROIs=nirs.util.convertlabels2roi(dataChannel(i).probe,obj.listOfROIs);
                else
                    listOfROIs=obj.listOfROIs;
                end
                if(isa(dataChannel(i),'nirs.core.Data'))
                   
                    dataROI(i,1)=nirs.util.roiAverage(dataChannel(i),listOfROIs,[],[],obj.weighted);
                else
                    if size(listOfROIs,1) == 1
                        [~,dataROI(i,1)]=nirs.util.roiAverage(dataChannel(i),listOfROIs,[],[],obj.weighted);
                    else
                        ROIS = {};
                        Names = {};
                        for ROI_idx = 1:size(listOfROIs,1)
                            source = listOfROIs(ROI_idx,:).source{1}';
                            detector = listOfROIs(ROI_idx,:).detector{1}';
                            ROIS{ROI_idx} = table(source, detector);
                            Names{ROI_idx} = listOfROIs(ROI_idx,:).Name{1};
                        end
                        [~,dataROI(i,1)]=nirs.util.roiAverage(dataChannel(i),ROIS,Names,[],obj.weighted);
                    end
                end
            end
            

        end
        
    end
end

