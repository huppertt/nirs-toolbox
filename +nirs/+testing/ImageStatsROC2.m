classdef ImageStatsROC2
%% ImageStatsROC - This class will perform ROC analysis. 
% Takes in a function to simulate data and an analysis pipeline that ends
% with a reconstructed ImageStats object.
% 
% Example 1:
%     % load some data
%     raw = nirs.io.loadDirectory( '~/resting_state', {} );
%     sd  = unique([raw(1).probe.link.source raw(1).probe.link.detector], 'rows');
%     n   = size(sd,1);
% 
%     % simulation function
%     rndData = @() raw(randi(length(raw)));
%     rndChan = @() sd(rand(n,1)<0.5,:);
%     simfunc = @() nirs.testing.simData( rndData(), [], 7, rndChan() );
% 
%     % setup ROC
%     test = nirs.testing.ChannelStatsROC();
%     test.simfunc = simfunc;
% 
%     % run ROC
%     test = test.run( 10 );
% 
%     % draw ROC curves
%     test.draw()
% 
% Example 2:
%     % pipeline
%     p = nirs.modules.Resample();
%     p = nirs.modules.OpticalDensity(p);
%     p = nirs.modules.BeerLambertLaw(p);
%     p = nirs.modules.AR_IRLS(p);
%     p.verbose = false;
%     p = nirs.modules.MixedEffects(p);
%     p.formula = 'beta ~ -1 + cond';
% 
%     test = nirs.testing.ChannelStatsROC(p, @nirs.testing.simDataSet);
% 
%     test = test.run(5);
% 
% Example 3:
%     % pipeline
%     p = nirs.modules.Resample();
%     p = nirs.modules.OpticalDensity(p);
%     p = nirs.modules.BeerLambertLaw(p);
%     p = nirs.modules.AR_IRLS(p);
%     p.verbose = false;
% 
%     p = nirs.modules.MixedEffects(p);
%     p.formula = 'beta ~ -1 + cond';
%     p.dummyCoding = 'full';
% 
%     % going to use ttest to subtract cond B from A
%     pipeline.run = @(data) p.run(data).ttest([1 -1]);
% 
%     % sim function
%     randStim = @(t) nirs.testing.randStimDesign(t, 2, 7, 2);
%     simfunc = @()nirs.testing.simDataSet([], [], randStim, [5 2]', []);
% 
%     % ROC test
%     test = nirs.testing.ChannelStatsROC(pipeline, simfunc);
% 
%     test = test.run(3);

    properties
        %simfunc  = @()nirs.testing.simDataImage2([], [], [], {'BA-46_R'}, [])
        simfunc  = @nirs.testing.simDataImage
        pipeline
    end
    
    properties (SetAccess = protected)
       truth
       XY
       XY_null
    end
    
    methods
        % constructor
        function obj = ImageStatsROC2( pipeline, simfunc )
           if nargin < 1
               p = nirs.modules.Resample();
               p = eeg.modules.BandPassFilter(p);
               p.lowpass = 0.5;
               p = nirs.modules.OpticalDensity(p);
               p = nirs.modules.AR_IRLS(p);
               p.verbose = false;
               p = nirs.modules.ImageReconMFX2(p);
               p.formula = 'beta ~ -1 + cond';
               
               
               obj.pipeline = p;
           else
               obj.pipeline = pipeline;
           end
           
           if nargin > 1
               obj.simfunc = simfunc;
           end
        end
        
        function obj = run(obj, iter)
            XY = [];
            XY_null = [];
            truth_all = [];
            for k = 1:iter
               [data, truth,fwdmodel,~,nulldata] = obj.simfunc();
               % pipeline stats
                for i=1:length(obj.pipeline)
                   
                   PL=nirs.modules.pipelineToList(obj.pipeline(i));
                   for j=1:length(PL)
                       if(isa(PL{j},'nirs.modules.ImageReconMFX2'))
                            PL{j}.mesh=fwdmodel.mesh;
                            PL{j}.probe('default')=fwdmodel.probe;
                            PL{j}.jacobian('default')=fwdmodel.jacobian('spectral');
                            %PL{j}.basis=nirs.inverse.basis.identity(length(fwdmodel.mesh.nodes));
                            PL{j}.basis=nirs.inverse.basis.gaussian(fwdmodel.mesh);
                            %PL{j}.basis=nirs.inverse.basis.freesurfer_wavelet(3);
                       end
                   end
                   obj.pipeline(i)=nirs.modules.listToPipeline(PL);
                   
                   XY = [XY, obj.pipeline(i).run(data)];
                   XY_null = [XY_null, obj.pipeline(i).run(nulldata)]; 
                   truth_all = [truth_all, truth];
                end 
                display(['Iter: ', num2str(k)]);
            end
            obj.XY = XY;
            obj.XY_null = XY_null;
            obj.truth = truth_all;
        end
    end    
end

