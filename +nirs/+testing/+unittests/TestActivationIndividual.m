classdef TestActivationIndividual < matlab.unittest.TestCase

    properties
        SubjStats
    end
    
    methods
        function obj = TestActivationIndividual()            
            
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
            
            link = link(randperm(height(link)),:);

            probe = nirs.core.Probe(srcPos, detPos, link);
            
            %% Create task design
            time = (0:2000) * (1/10);
            
            stim = Dictionary();
            
            stim('Low') = nirs.design.StimulusEvents('Low',[25 97 149],[17 6 12],[1 1 1]);
            stim('Medium') = nirs.design.StimulusEvents('Medium',[5 80 112],[12 15 16],[1 1 1]);
            stim('High') = nirs.design.StimulusEvents('High',[50 60 134],[6 10 7],[1 1 1]);
            
            beta = [20 60 100]';
            
            %% Generate data
            rng('default');
            noise = nirs.testing.simARNoise(probe,time);
            raw = nirs.testing.simData(noise,stim,beta);
            
            %% Preprocess
            job_preproc = nirs.modules.OpticalDensity();
            job_preproc = nirs.modules.Resample(job_preproc);
            job_preproc.Fs = 4.253;
            job_preproc = nirs.modules.BeerLambertLaw(job_preproc);
            
            hb = job_preproc.run(raw);
            
            %% Calculate activation
            canon = nirs.design.basis.Canonical();
            canon.incDeriv = true;
            
            basis = Dictionary();
            basis('default') = canon;
            
            job_act = nirs.modules.GLM();
            job_act.basis = basis;
            job_act.trend_func = @(t) nirs.design.trend.dctmtx(t,1/128);
            
            obj.SubjStats = job_act.run(hb);
            
        end
    end
    methods (Test)
        
        function testConditionHighGreaterThanMedium( obj )
            tbl = obj.SubjStats.ttest('High:01-Medium:01').table;
            tstat_active_hbo = tbl.tstat( tbl.source==1 & strcmp(tbl.type,'hbo'));
            tstat_inactive_hbo = tbl.tstat( tbl.source==2 & strcmp(tbl.type,'hbo'));
            tstat_active_hbr = tbl.tstat( tbl.source==1 & strcmp(tbl.type,'hbr'));
            tstat_inactive_hbr = tbl.tstat( tbl.source==2 & strcmp(tbl.type,'hbr'));
            obj.verifyGreaterThan( min(tstat_active_hbo) , 3 );
            obj.verifyGreaterThan( min(-tstat_active_hbr) , 3 );
            obj.verifyLessThan( max(abs(tstat_inactive_hbo)) , 3 );
            obj.verifyLessThan( max(abs(tstat_inactive_hbr)) , 3 );    
        end
        
        function testConditionMediumGreaterThanLow( obj )
            tbl = obj.SubjStats.ttest('Medium:01-Low:01').table;
            tstat_active_hbo = tbl.tstat( tbl.source==1 & strcmp(tbl.type,'hbo'));
            tstat_inactive_hbo = tbl.tstat( tbl.source==2 & strcmp(tbl.type,'hbo'));
            tstat_active_hbr = tbl.tstat( tbl.source==1 & strcmp(tbl.type,'hbr'));
            tstat_inactive_hbr = tbl.tstat( tbl.source==2 & strcmp(tbl.type,'hbr'));
            obj.verifyGreaterThan( min(tstat_active_hbo) , 3 );
            obj.verifyGreaterThan( min(-tstat_active_hbr) , 3 );
            obj.verifyLessThan( max(abs(tstat_inactive_hbo)) , 3 );
            obj.verifyLessThan( max(abs(tstat_inactive_hbr)) , 3 );
        end

        function testConditionHighGreaterThanLow( obj )
            tbl = obj.SubjStats.ttest('High:01-Low:01').table;
            tstat_active_hbo = tbl.tstat( tbl.source==1 & strcmp(tbl.type,'hbo'));
            tstat_inactive_hbo = tbl.tstat( tbl.source==2 & strcmp(tbl.type,'hbo'));
            tstat_active_hbr = tbl.tstat( tbl.source==1 & strcmp(tbl.type,'hbr'));
            tstat_inactive_hbr = tbl.tstat( tbl.source==2 & strcmp(tbl.type,'hbr'));
            obj.verifyGreaterThan( min(tstat_active_hbo) , 3 );
            obj.verifyGreaterThan( min(-tstat_active_hbr) , 3 );
            obj.verifyLessThan( max(abs(tstat_inactive_hbo)) , 3 );
            obj.verifyLessThan( max(abs(tstat_inactive_hbr)) , 3 );
        end
        
    end
end