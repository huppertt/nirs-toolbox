classdef NIRS_SPM_GLM < nirs.modules.AbstractGLM
    %% This module is a wrapper for the NIRS-SPM code
    %
    % Options:
    %     basis       - a Dictionary object containing temporal bases using stim name as key
    %     verbose     - flag to display progress
    %     trend_func  - a function that takes in a time vector and returns trend regressors
    
    %     spm_hpf -  {<wavelet> or <DCT,##> }
    %     spm_lpf -  {<none> or <hrf>  or <gaussian,FWHM>  }
    %     spm_cVi -  {<none> or <AR(##)> }
    
    properties
        spm_hpf='DCT,128';
        spm_lpf='none';
        spm_cVi='AR(1)';
    end
    
    methods
        function obj = NIRS_SPM_GLM( prevJob )
            if nargin > 0, obj.prevJob = prevJob; end
            
            obj.name = 'GLM via NIRS-SPM';
            obj.basis('default') = nirs.design.basis.Canonical();
            
            if(isempty(which('spm_filter_HPF_LPF_WMDL')))
                warning('NIRS-SPM not found on Matlab Path');
            end
            obj.citation='Ye, J. C., Tak, S., Jang, K. E., Jung, J., & Jang, J. (2009). NIRS-SPM: statistical parametric mapping for near-infrared spectroscopy. Neuroimage, 44(2), 428-447.';
            
            
        end
        
        function S = runThis( obj, data )
            vec = @(x) x(:);
            
            if(isempty(which('spm_filter_HPF_LPF_WMDL')))
                error('NIRS-SPM not found on Matlab Path');
            end
            
            for i = 1:length(data)
                % get data
                d  = data(i).data;
                t  = data(i).time;
                Fs = data(i).Fs;
                
                probe = data(i).probe;
                
                % make sure data is in order
                [probe.link, idx] = sortrows(probe.link, {'source', 'detector','type'});
                d = d(:, idx);
                
                % get experiment design
                [X, names] = obj.createX( data(i) );
                C = obj.getTrendMatrix( t );
                
                for id=1:size(C,2)
                    cnames{id}=sprintf('trend_%d',i);
                end
                
                % check model
                obj.checkRank( [X C] )
                obj.checkCondition( [X C] )
                
                if(rank([X C]) < size([X C],2) & obj.goforit)
                    disp('Using PCA regression model');
                    [U,s,V]=nirs.math.mysvd([X C]);
                    lst=find(diag(s)>eps(1)*10);
                    V=V(:,lst);
                    SPM.xX.X      = U(:,lst)*s(lst,lst);
                    
                    
                else
                    SPM.xX.X      = [X C];
                end
                
                nscan=length(t);
                
                % Create the SPM structure
                
                SPM.xX.iH     = [];
                SPM.xX.iC     = [1:size(X,2)];
                SPM.xX.iB     = [1:size(C,2)] + size(X,2);
                SPM.xX.iG     = [];
                SPM.xX.name   = {names{:} cnames{:}};
                
                str='Empty dummy file';
                save('temp.mat','str')
                SPM.nirs.fname='temp.mat';
                
                % High-pass filtering
                if isempty(strfind(obj.spm_hpf, 'wavelet')) == 0 % wavelet-MDL
                    SPM.xX.K.HParam.type = 'Wavelet-MDL';
                elseif isempty(strfind(obj.spm_hpf, 'DCT')) == 0 % DCT
                    index_cutoff = find(obj.spm_hpf == ',');
                    if isempty(index_cutoff) == 1
                        cutoff = 128;
                    else
                        cutoff = str2num(obj.spm_hpf(index_cutoff+1:end));
                    end
                    SPM.xX.K.HParam.type = 'DCT';
                    SPM.xX.K.HParam.M = cutoff;
                else
                    SPM.xX.K.HParam.type = 'none';
                end
                
                
                % Low-pass filtering
                
                if isempty(strfind(obj.spm_lpf, 'hrf')) == 0 % hrf smoothing
                    SPM.xX.K.LParam.type = 'hrf';
                elseif isempty(strfind(obj.spm_lpf, 'gaussian')) == 0 % Gaussian smoothing
                    index_FWHM = find(obj.spm_lpf == ',');
                    if isempty(index_FWHM) == 1
                        FWHM = 4;
                    else
                        FWHM = str2num(obj.spm_lpf(index_FWHM+1:end));
                    end
                    SPM.xX.K.LParam.FWHM = FWHM;
                    SPM.xX.K.LParam.type = 'Gaussian';
                else
                    SPM.xX.K.LParam.type = 'none';
                end
                
                
                
                K = struct( 'HParam', SPM.xX.K.HParam,...
                    'row', [1:nscan]',...
                    'RT', 1/Fs,...
                    'LParam', SPM.xX.K.LParam);
                SPM.xX.K = spm_filter_HPF_LPF_WMDL(K);
                
                if(strcmp(obj.spm_cVi,'none'))
                    SPM.xVi.V  = speye(sum(nscan));
                    cVi        = 'i.i.d';
                else
                    SPM.xVi.Vi = spm_Ce(nscan,0.2);
                    cVi        = 'AR(0.2)';
                end
                
                SPM.xVi.form = cVi;
                
                data(i)=data(i).sorted({'type','source','detector'});
                probe=data(i).probe;
                utype=unique(data(i).probe.link.type);
                covb=[];
                beta=[];
                SPMorig=SPM;
                for tyIdx=1:length(utype)
                    lst=ismember(data(i).probe.link.type,utype{tyIdx});
                    Y=data(i).data(:,lst);
                    SPM=SPMorig;
                    if isfield(SPM.xVi, 'V') == 1 % precoloring method
                        SPM = rmfield(SPM, 'xVi');
                        [SPM] = precoloring(SPM, Y);
                    elseif isfield(SPM.xVi, 'V') == 0
                        [SPM] = prewhitening(SPM, Y, pwd);
                    end
                    beta = [beta; vec(SPM.nirs.beta(SPM.xX.iC,:))];
                    covb = blkdiag(covb,kron(SPM.xX.Bcov(SPM.xX.iC,SPM.xX.iC),SPM.nirs.ResSS./SPM.xX.trRV));
                end
                
                %
                % put stats
                ncond = length(names);
                nchan = size(data(i).probe.link, 1);
                
                link = repmat( probe.link, [ncond 1] );
                cond = repmat(names(:)', [nchan 1]);
                cond = cond(:);
                
                if(isempty(~strfind(class(probe),'nirs')))
                    S(i) = nirs.core.ChannelStats();
                elseif(isempty(~strfind(class(probe),'eeg')))
                    S(i) = eeg.core.ChannelStats();
                else
                    warning('unsupported data type');
                    S(i) = nirs.core.ChannelStats();
                end
                
                S(i).variables = [link table(cond)];
                S(i).beta = beta;
                
                
                S(i).covb = covb;
                
                if(rank([X C]) < size([X C],2) & obj.goforit)
                    stats.beta=V*stats.beta;
                    for j=1:size(stats.covb,3)
                        c(:,:,j)=V*squeeze(stats.covb(:,:,j))*V';
                    end
                    stats.covb=c;
                end
                
                S(i).dfe  = SPM.xX.erdf;
                
                S(i).description = data(i).description;
                
                S(i).demographics   = data(i).demographics;
                S(i).probe          = probe;
                
                stim=Dictionary;
                for j=1:data(i).stimulus.count;
                    ss=data(i).stimulus.values{j};
                    if(isa(ss,'nirs.design.StimulusEvents'))
                        s=nirs.design.StimulusEvents;
                        s.name=ss.name;
                        s.dur=mean(ss.dur);
                        stim(data(i).stimulus.keys{j})=s;
                    end
                end
                
                S(i).basis.base=obj.basis;
                S(i).basis.Fs=Fs;
                S(i).basis.stim=stim;
                
                
                % print progress
                obj.printProgress( i, length(data) )
                
            end
            delete('temp.mat');
            
        end
        
    end
end


