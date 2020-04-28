classdef ImageReconMFX2 < nirs.modules.AbstractModule
    %This is the mixed effects image reconstruction model
    %   This model preforms single-subject or group-level image
    %   reconstruction using ReML
    
    properties
        formula = 'beta ~ -1 + cond*group + (1|subject)';
        jacobian = Dictionary(); % key is subject name or "default"
        dummyCoding = 'full';
        centerVars=false
        
        basis;  % Basis set from nirs.inverse.basis
        mask = [];   % Mask (e.g. cortical contraint)
        prior = Dictionary();  % The prior on the image; default = 0 (Min Norm Estimate)
        mesh;  % Reconstruction mesh (needed to create the ImageStats Class)
        probe = Dictionary();  % This is the probe for the Jacobian model.  These needs to match the jacobian      
    end
    
    methods
        
        function obj = ImageReconMFX2( prevJob )
            obj.name = 'Image Recon w/ Random Effects';
            if nargin > 0
                obj.prevJob = prevJob;
            end
                      
            obj.citation{1}='Abdelnour, F., B. Schmidt, and T. J. Huppert. "Topographic localization of brain activation in diffuse optical imaging using spherical wavelets." Physics in medicine and biology 54.20 (2009): 6383.';
            obj.citation{2}='Abdelnour, F., & Huppert, T. (2011). A random-effects model for group-level analysis of diffuse optical brain imaging. Biomedical optics express, 2(1), 1-25.';
            obj.citation{3}='Abdelnour, F., Genovese, C., & Huppert, T. (2010). Hierarchical Bayesian regularization of reconstructions for diffuse optical tomography using multiple priors. Biomedical optics express, 1(4), 1084-1103.';
            
        end
        
        function G = runThis( obj,S )
            
            %rescale=(50*length(unique(nirs.getStimNames(S)))*length(S));
            
            % demographics info
            demo = nirs.createDemographicsTable( S );
            if(isempty(demo))
                for idx=1:length(S);
                    S(idx).demographics('subject')='default';
                end
            end
            if(~ismember(demo.Properties.VariableNames,'subject'));
                for idx=1:length(S);
                    S(idx).demographics('subject')='default';
                end
            end
            
            
            % center numeric variables
            if obj.centerVars
                n = demo.Properties.VariableNames;
                for i = 1:length(n)
                    if all( isnumeric( demo.(n{i}) ) )
                        demo.(n{i}) = demo.(n{i}) - mean( demo.(n{i}) );
                    end
                end
            end
            
            
            % Mask of the reconstruction volume
            if(~isempty(obj.mask))
               LstInMask=find(obj.mask);
            else
               LstInMask=1:size(obj.basis.fwd,1);
            end
                        
            %% Wavelet ( or other tranform )
            Basis = obj.basis.fwd; % W = W(1:2562,:);
            
            %Let's make the forward models
            L = obj.jacobian;
            Lfwdmodels=Dictionary();
            Probes=Dictionary();
            
            for i = 1:L.count
                key = L.keys{i};
                J =L(key);
                flds=fields(J);
                
                if(~ismember(key,obj.probe.keys))
                     Probes(key)=obj.probe('default');
                else
                     Probes(key)=obj.probe(key);
                end
                l=[];
                for j=1:length(flds)
                    l=setfield(l,flds{j},(J.(flds{j})(:,LstInMask).*(J.(flds{j})(:,LstInMask)>10*eps(1)))*...
                        Basis(LstInMask,:));
                end
                Lfwdmodels(key)=l;
            end   
            
            %Make sure the probe and data link match
            for idx=1:length(S)
                sname=S(idx).demographics('subject');
                if(ismember(sname,Lfwdmodels.keys))
                    key=sname;
                else
                    key='default';
                end
                thisprobe=Probes(key);
                [ia,ib]=ismember(thisprobe.link,S(idx).probe.link,'rows');
                ibAll=[];
                for i=1:length(S(idx).conditions);
                    ibAll=[ibAll; ib+(length(ib)*(i-1))];
                end   
                S(idx).variables=S(idx).variables(ibAll,:);
                S(idx).beta=S(idx).beta(ibAll);
                S(idx).covb=S(idx).covb(ibAll,ibAll);
                S(idx).probe.link=S(idx).probe.link(ib,:);
                
            end            
            G.l = J;
            G.beta = S.table;
            G.covb = S.covb;
        end
        
    end
    
end

