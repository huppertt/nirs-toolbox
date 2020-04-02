classdef ApplyROIMFX < nirs.modules.AbstractModule
%% ApplyROIMFX - Performs ROI averaging by creating a within-subject MFX 
%                model with channel ID (or source and detector IDs) treated as random effects
% Usage:
% job = advanced.nirs.modules.ApplyROIMFX();
% job.listOfROIs(1,:) = array2table({[1 2 2 3 4 5 5 6],[1 1 2 2 3 3 4 4],'PFC'});
% dataROI = job.run( dataChannel );
%
% Options: 
%     listOfROIs - [#ROI x 3] table of ROIs (source, detector, name)
%     model_optodes - flag to enable modeling source and detector IDs rather than channel IDs [default: false]
%     covariance_pattern - See fitlme() documentation. [default: 'Isotropic']
%
% See also advanced.nirs.util.threshold2ROI() for creating ROI tables from significance masks

    properties
        listOfROIs = table({},{},{},'VariableNames',{'source','detector','name'});
        model_optodes = false;
        covariance_pattern = 'Isotropic';
    end
    
    methods
        function obj = ApplyROIMFX( prevJob )
         	if nargin > 0, obj.prevJob = prevJob; end
            obj.name = 'Apply ROI MFX';
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

            if obj.model_optodes
                formula = 'beta ~ -1 + cond + (1|source) + (1|detector)';
                covpattern = {obj.covariance_pattern, obj.covariance_pattern};
            else
                formula = 'beta ~ -1 + cond + (1|channel)';
                covpattern = obj.covariance_pattern;
            end

            dataChannel = dataChannel.sorted();
            nROI = height(probe.link);
            
            dataROI = dataChannel;
            switch class(dataChannel)
                    
                case {'nirs.core.ChannelStats'}
                    
                    vars = dataChannel(1).variables;
                    for i = 1:height(vars)
                        channel{i,1} = sprintf('s%id%i',vars.source(i),vars.detector(i));
                    end
                    datatbl = [vars table(dataChannel(1).beta,'VariableNames',{'beta'})];
                    datatbl.source = categorical(datatbl.source);
                    datatbl.detector = categorical(datatbl.detector);
                    datatbl.channel = channel;
                    
                    types = unique(vars.type,'stable');
                    ncond = length(dataChannel(1).conditions);
                    if isnumeric(types)
                        newvars = table({},{},[],{},{},'VariableNames',{'source','detector','type','ROI','cond'});
                    else
                        newvars = table({},{},{},{},{},'VariableNames',{'source','detector','type','ROI','cond'});
                    end
                    
                    for i = 1:length(dataROI)
                        dataROI(i).probe = probe;
                        dataROI(i).variables = newvars;
                        dataROI(i).beta = zeros(nROI*ncond,1);
                        dataROI(i).covb = zeros(nROI*ncond);
                    end

                    msg = 0;
                    for i = 1:nROI
                        
                        fprintf(repmat('\b',1,msg));
                        msg = fprintf('Processing ROI %i/%i ',i,nROI);
                        
                        sources = probe.link.source{i};
                        detectors = probe.link.detector{i};
                        type = probe.link.type(i);
                            
                        % set output index
                        out_idx = (i-1)*ncond + (1:ncond);

                        % find relevant channel indices for ROI/type
                        incl_chan = false(height(vars),1);
                        if isnumeric(type)
                            has_t = (vars.type == type);
                        else
                            has_t = strcmpi(vars.type,type);
                        end
                        for k = 1:length(sources)
                            has_s = (vars.source == sources(k));
                            has_d = (vars.detector == detectors(k));
                            incl_chan = incl_chan | ( has_s & has_d & has_t );
                        end

                        % Create MFX model
                        mdl0 = fitlme(datatbl,formula,'DummyVarCoding','full','Exclude',~incl_chan);

                        % Get design matrices
                        X = mdl0.designMatrix('fixed');
                        Z = mdl0.designMatrix('random');
                        cond = strrep(mdl0.CoefficientNames(:),'cond_','');
                        tmpvar = [repmat(probe.link(i,:),ncond,1) table(cond)];

                        % Run model for each subject
                        msg2 = 0;
                        for sub = 1:length(dataChannel)

                            assert(isequal(dataChannel(1).variables,dataChannel(sub).variables),'Variable table mismatch');

                            fprintf(repmat('\b',1,msg2));
                            msg2 = fprintf('(%4.4f%%)',100*(sub-1)/length(dataChannel));
                            
                            try
                                mdl = fitlmematrix(X,dataChannel(sub).beta,Z,[],'CovariancePattern',covpattern,'FitMethod','ML');

                                dataROI(sub).variables(out_idx,:) = tmpvar;
                                dataROI(sub).beta(out_idx) = mdl.Coefficients.Estimate;
                                dataROI(sub).covb(out_idx,out_idx) = mdl.CoefficientCovariance;
                            catch err
                                warning(err.message);
                            end

                        end
                        fprintf(repmat('\b',1,msg2));
                    end
                    fprintf(repmat('\b',1,msg));
                    
                otherwise
                    error('Type %s not implemented.',class(dataChannel));
            end
        end            
        
        % Generates a new probe for the ROIs (source & detector fields in
        % link are now arrays and link has 'ROI' column with the region name)
        function probeROI = get_probeROI( obj , probe )
            
            source = obj.listOfROIs.source;
            detector = obj.listOfROIs.detector;
            name = obj.listOfROIs.name;
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

    end
end

