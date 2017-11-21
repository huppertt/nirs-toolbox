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
%     weighted   - logical (False); Option to enable variance-weighted averaging
%                                   (currently only supported for time series)

    properties
        listOfROIs = table({},{},{},'VariableNames',{'sources','detectors','names'});
        weighted = false;
    end
    
    methods
        function obj = ApplyROI( prevJob )
         	if nargin > 0, obj.prevJob = prevJob; end
            obj.name = 'Apply ROI';
        end
        
        % Apply ROI averaging to data and update probe
        function dataROI = runThis( obj , dataChannel )
            
            oldprobe = dataChannel(1).probe;
            for i = 1:length(dataChannel)
                if ~isequal(oldprobe,dataChannel(i).probe)
                    error('Please standardize probe between scans');
                end
            end
            
            probe = obj.get_probeROI(oldprobe);
            
            if isempty(probe)
                error('Problem setting up channel probe and ROIs');
            end
            if ~any(strcmp(probe.link.Properties.VariableNames,'ROI'))
                error('No ROIs detected in probeROI');
            end
            
            if obj.weighted && ~any(strcmp({'nirs.core.Data','nirs.core.ChannelStats'},class(dataChannel)))
                error('Weighted option is only supported for these classes: nirs.core.Data, nirs.core.ChannelStats');
            end
            
            dataROI = dataChannel;
            switch class(dataChannel)
                case {'nirs.core.Data'}
                    projmat = obj.getMapping(oldprobe,probe,true);
                    for i = 1:length(dataChannel)
                        tmpprojmat = projmat;
                        if obj.weighted
                            % Weights are the inverse of the robust standard deviation.  
                            % We apply to the mapping matrix and rescale so each column/ROI sums to 1
                            weights = 1 ./ (1.4826 * mad( dataChannel(i).data ));
                            tmpprojmat = bsxfun( @times , tmpprojmat , weights' );
                            tmpprojmat = bsxfun( @rdivide , tmpprojmat , sum(tmpprojmat) );
                        end
                        dataROI(i).probe = probe;
                        dataROI(i).data = dataChannel(i).data * tmpprojmat;
                    end
                    
                case {'nirs.core.ChannelStats'}
                    projmat = obj.getMapping(oldprobe,probe,true);
                    for i = 1:length(dataChannel)
                        dataROI(i).probe = probe;
                        vars = dataChannel(i).variables;
                        [vars,sort_inds] = sortrows( vars , {'cond','source','detector','type'} );
                        conds = unique(vars.cond,'stable');
                        numcond = length(conds);
                        condvec = repmat(conds(:)',[height(probe.link) 1]);
                        condprojmat = kron(eye(numcond),projmat);
                        
                        beta = dataChannel(i).beta(sort_inds);
                        covb = dataChannel(i).covb(sort_inds,sort_inds);
                        
                        if obj.weighted
                            [u, s, ~] = svd(covb, 'econ');
                            w = u * pinv(s).^.5 * u';
                            condprojmat = w * condprojmat;
                            condprojmat = bsxfun(@rdivide,condprojmat,sum(condprojmat));
                        end
                        
                        dataROI(i).beta = condprojmat' * beta;
                        dataROI(i).covb = condprojmat' * covb * condprojmat;
                        dataROI(i).variables = repmat(probe.link,numcond,1);
                        dataROI(i).variables.cond = condvec(:);
                    end
                    
                case {'nirs.core.sFCStats'}
                    projmat = obj.getMapping(oldprobe,probe,false);
                    for i = 1:length(dataChannel)
                        
                        dataROI(i).probe = probe;
                        numconds = length(dataChannel(i).conditions);
                        ROI_size = [size(projmat,2) size(projmat,2) numconds];
                        goodvals = ~isnan(dataChannel(i).R) & isfinite(dataChannel(i).Z);
                        dataChannel(i).R(~goodvals) = 0;
                        dataROI(i).R = zeros(ROI_size);
                        
                        % Average Z-transformed R-values
                        for j = 1:length(dataChannel(i).conditions)
                            Z = projmat' * dataChannel(i).Z(:,:,j) * projmat; % ROI sum of channel Z-values
                            numnotnan = projmat' * goodvals(:,:,j) * projmat; % ROI sum of channel NaNs
                            Z = Z ./ numnotnan; % Sum to mean of non-nanvals
                            dataROI(i).R(:,:,j) = tanh( Z ); % Z-to-R
                        end
                        
                        % Average StdErr
                        if ~isempty(dataChannel(i).ZstdErr)
                            
                            dataChannel(i).ZstdErr(~goodvals) = 0;
                            dataROI(i).ZstdErr = zeros(ROI_size);
                                                   
                            for j = 1:length(dataChannel(i).conditions)
                                ZstdErr = projmat' * dataChannel(i).ZstdErr(:,:,j).^2 * projmat; % ROI sum of channel values
                                numnotnan = projmat' * goodvals(:,:,j) * projmat; % ROI sum of channel NaNs
                                dataROI(i).ZstdErr(:,:,j) = sqrt(ZstdErr ./ numnotnan); % Mean of non-nanvals
                            end
                            
                        end
                        
                        % Maintain symmetry
                        dataROI(i).R = (dataROI(i).R + permute(dataROI(i).R,[2 1 3])) ./ 2;
                        dataROI(i).ZstdErr = (dataROI(i).ZstdErr + permute(dataROI(i).ZstdErr,[2 1 3])) ./ 2;
                        
                    end
                    
                case {'nirs.core.sFCBetaStats'}
                    projmat = obj.getMapping(oldprobe,probe,false);
                    for i = 1:length(dataChannel)
                        
                        dataROI(i).probe = probe;
                        numconds = length(dataChannel(i).conditions);
                        ROI_size = [size(projmat,2) size(projmat,2) numconds];
                        goodvals = ~isnan(dataChannel(i).beta);
                        dataChannel(i).beta(~goodvals) = 0;
                        dataROI(i).beta = zeros(ROI_size);
                        
                        % Average Z-transformed R-values
                        for j = 1:length(dataChannel(i).conditions)
                            beta = projmat' * dataChannel(i).beta(:,:,j) * projmat; % ROI sum of channel Z-values
                            numnotnan = projmat' * goodvals(:,:,j) * projmat; % ROI sum of channel NaNs
                            beta = beta ./ numnotnan; % Sum to mean of non-nanvals
                            dataROI(i).beta(:,:,j) = beta; % Z-to-R
                        end
                        
                        % Average StdErr
                        if ~isempty(dataChannel(i).covb)
                            
                            goodvals = ~isnan(dataChannel(i).covb);
                            dataChannel(i).covb(~goodvals) = 0;
                            dataROI(i).covb = zeros(ROI_size);
                                                   
                            for j = 1:length(dataChannel(i).conditions)
                                for k = 1:length(dataChannel(i).conditions)
                                    covb = projmat' * dataChannel(i).covb(:,:,j,k) * projmat; % ROI sum of channel values
                                    numnotnan = projmat' * goodvals(:,:,j,k) * projmat; % ROI sum of channel NaNs
                                    dataROI(i).covb(:,:,j,k) = covb ./ numnotnan; % Mean of non-nanvals
                                end
                            end
                            
                        end
                        
                        % Maintain symmetry
                        dataROI(i).beta = (dataROI(i).beta + permute(dataROI(i).beta,[2 1 3])) ./ 2;
                        dataROI(i).covb = (dataROI(i).covb + permute(dataROI(i).covb,[2 1 3])) ./ 2;
                        
                    end
                    
                otherwise
                    error('Type %s not implemented.',class(dataChannel));
            end
        end            
        
        % Generates a new probe for the ROIs (source & detector fields in
        % link are now arrays and link has 'ROI' column with the region name)
        function probeROI = get_probeROI( obj , probe )
            
            source = obj.listOfROIs.sources;
            detector = obj.listOfROIs.detectors;
            name = obj.listOfROIs.names;
            link = probe.link;
            types = unique(link.type,'stable');
            if isnumeric(types)
                link = table({},{},[],{},'VariableNames',{'source','detector','type','ROI'});
            else
                link = table({},{},{},{},'VariableNames',{'source','detector','type','ROI'});
            end
            for i = 1:length(source)
                for j = 1:length(types)
                    link = [link; table(source(i),detector(i),types(j),name(i),'VariableNames',{'source','detector','type','ROI'})];
                end
            end

            if any(strcmp(probe.link.Properties.VariableNames,'hyperscan'))
                inds_A = strfind( probe.link.hyperscan' , 'A' );
                inds_B = strfind( probe.link.hyperscan' , 'B' );
                hyper_source_offset = min(probe.link.source(inds_B)) - min(probe.link.source(inds_A));
                hyper_detector_offset = min(probe.link.detector(inds_B)) - min(probe.link.detector(inds_A));
                linkA = link; linkA.hyperscan = repmat('A',[height(link) 1]);
                linkB = link; linkB.hyperscan = repmat('B',[height(link) 1]);
                for i = 1:height(linkB)
                    linkB.source{i} = linkB.source{i} + hyper_source_offset;
                    linkB.detector{i} = linkB.detector{i} + hyper_detector_offset;
                end
                link = [linkA; linkB];
            end
            
            probe.link = link;
            probeROI = probe;
        end
                
        % Returns a [#channel x 1] logical of channels within specified ROI
        function inds = getChannelInds(obj,probeChannel,s,d,t)
            chanlink = probeChannel.link;
            inds = false(height(chanlink),1);
            for i = 1:length(s)
                if isnumeric(t)
                    inds = inds | (chanlink.source==s(i) & chanlink.detector==d(i) & chanlink.type == t);
                else
                    inds = inds | (chanlink.source==s(i) & chanlink.detector==d(i) & strcmp(chanlink.type,t));
                end
            end
        end
        
        % Returns a [#channel x #ROI] binary projection matrix
        function mapping = getMapping( obj , probeChannel , probeROI , scale )
            if nargin<2, scale = false; end
            chanlink = probeChannel.link;
            roilink = probeROI.link;
            num_chan = height(chanlink);
            num_ROI = height(roilink);
            
            mapping = zeros(num_chan,num_ROI);
            for i = 1:num_ROI
                inds = obj.getChannelInds(probeChannel,roilink.source{i},roilink.detector{i},roilink.type(i));
                if scale
                    mapping(inds,i) = 1/sum(inds);
                else
                    mapping(inds,i) = 1;
                end
            end
        end
        
    end
end

