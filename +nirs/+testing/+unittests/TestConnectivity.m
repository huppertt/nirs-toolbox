classdef TestConnectivity < matlab.unittest.TestCase
    
    properties
        GroupStats
    end
    
    methods
        function obj = TestConnectivity()
            %% Create probe
            srcPos = [0 0 0;
                30 30 0];
            
            detPos = [ 0 30 0;
                30 0 0];
            
            link = [1 1 690;
                1 1 830;
                1 2 690;
                1 2 830;
                2 1 690;
                2 1 830;
                2 2 690;
                2 2 830 ];
            
            link = table(link(:,1), link(:,2), link(:,3), ...
                'VariableNames', {'source', 'detector', 'type'});
            
            probe = nirs.core.Probe(srcPos, detPos, link);
            
            %% Create task design
            time = (0:2000) * (1/10);
            
            stim = Dictionary();
            
            stim('Low') = nirs.design.StimulusEvents('Low',[25 97 149],[17 6 12],[1 1 1]);
            stim('Medium') = nirs.design.StimulusEvents('Medium',[5 80 112],[12 15 16],[1 1 1]);
            stim('High') = nirs.design.StimulusEvents('High',[50 60 134],[6 10 7],[1 1 1]);
            
            %% Generate data
            rng('default');
            for i = 1:4
                
                raw(i) = nirs.testing.simARNoise(probe,time);
                raw(i).stimulus = stim;
                raw(i).demographics('subject') = num2str(i);
                [X,names] = raw(i).getStimMatrix;
                
                
                % whiten channels
                means = mean(raw(i).data);
                raw(i).data = bsxfun( @minus , raw(i).data , means );
                C = cov(raw(i).data);
                [E,D] = eig( C );
                proj =  inv(sqrt(D)) * E';
                raw(i).data = (proj * raw(i).data')';

                % add correlations (positive load-related correlation
                % between S1-D1 and S2-D2)
                for j = 1:3
                    inds = find(X(:,j));
                    data_poscorr_wl1 = mean(raw(i).data(inds,[1 7]),2);
                    data_poscorr_wl2 = mean(raw(i).data(inds,[2 8]),2);
                    switch names{j}
                        case 'Low'
                            corr_weight = .1;
                        case 'Medium'
                            corr_weight = .5;
                        case 'High'
                            corr_weight = .9;
                    end
                    raw(i).data(inds,1) = corr_weight * data_poscorr_wl1 + (1-corr_weight) * raw(i).data(inds,1);
                    raw(i).data(inds,2) = corr_weight * data_poscorr_wl2 + (1-corr_weight) * raw(i).data(inds,2);
                    raw(i).data(inds,7) = corr_weight * data_poscorr_wl1 + (1-corr_weight) * raw(i).data(inds,7);
                    raw(i).data(inds,8) = corr_weight * data_poscorr_wl2 + (1-corr_weight) * raw(i).data(inds,8);                    
                end
                raw(i).data = bsxfun( @plus , raw(i).data , means );
            end
                        
            %% Randomize probe link sequence
            rand_inds = randperm(height(raw(1).probe.link));
            for i = 1:length(raw)
                raw(i).probe.link = raw(i).probe.link(rand_inds,:);
                raw(i).data = raw(i).data(:,rand_inds);
            end            
            
            %% Preprocess
            job = nirs.modules.OpticalDensity();
            job = nirs.modules.Resample(job);
            job.Fs = 4;
            job = nirs.modules.BeerLambertLaw(job);
            
            hb = job.run(raw);
            
            %% Level-1 connectivity
            job = nirs.modules.Connectivity();
            job.ignore = 0;
            job.min_event_duration = 1;
            job.divide_events = 1;
            SubjStats = job.run(hb);
            
            %% Level-2 connectivity
            job = nirs.modules.MixedEffectsConnectivity();
            job.formula = 'beta ~ -1 + cond + (1|subject)';
            job.robust = true;
            obj.GroupStats = job.run( SubjStats );
        end
    end
    
    methods (Test)
        
        function testConnectivityStats( obj )
            
            results = obj.GroupStats.table;
            
            active_inds = ( results.SourceOrigin==1 & results.DetectorOrigin==1 ...
                & results.SourceDest==2 & results.DetectorDest==2 ...
                & strcmp(results.TypeOrigin,results.TypeDest) );
            
            inactive_inds = ( results.SourceOrigin==1 & results.DetectorOrigin==2 ...
                & results.SourceDest==2 & results.DetectorDest==1 ...
                & strcmp(results.TypeOrigin,results.TypeDest) );
            
            tstat_active_low = results.t(active_inds & strcmp(results.condition,'Low'));
            tstat_active_medium = results.t(active_inds & strcmp(results.condition,'Medium'));
            tstat_active_high = results.t(active_inds & strcmp(results.condition,'High'));
            
            tstat_inactive_low = results.t(inactive_inds & strcmp(results.condition,'Low'));
            tstat_inactive_medium = results.t(inactive_inds & strcmp(results.condition,'Medium'));
            tstat_inactive_high = results.t(inactive_inds & strcmp(results.condition,'High'));
            
            obj.verifyLessThan( max(abs(tstat_inactive_low)) , 4 );
            obj.verifyLessThan( max(abs(tstat_inactive_medium)) , 4 );
            obj.verifyLessThan( max(abs(tstat_inactive_high)) , 4 );
            
            obj.verifyLessThan( max(abs(tstat_active_low)) , 4 );
            
            obj.verifyGreaterThan( min(tstat_active_medium) , 4 );
            obj.verifyLessThan( max(tstat_active_medium) , 16 );
            
            obj.verifyGreaterThan( min(tstat_active_high) , 16 );

        end
        
    end
end
            
            