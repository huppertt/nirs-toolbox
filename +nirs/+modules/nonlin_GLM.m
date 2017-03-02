classdef nonlin_GLM < nirs.modules.AbstractGLM
    %% AR_IRLS - Performs first-level per file GLM analysis.
    %
    % Options:
    %     basis       - a Dictionary object containing temporal bases using stim name as key
    %     verbose     - flag to display progress
    %     trend_func  - a function that takes in a time vector and returns trend regressors
    %
    % Example:
    %     j = nirs.modules.nonlin_GLM();
    %
    %     b = Dictionary();
    %     b('default') = nirs.design.basis.Canonical(); % default basis
    %     b('A')       = nirs.design.basis.Gamma();     % a different basis for condition 'A'
    %
    %     j.basis = b;
    %
    %     j.trend_func = @(t) nirs.design.trend.legendre(t, 3); % 3rd order polynomial for trend
    %
    % Note:
    %     trend_func must at least return a constant term unless all baseline periods are
    %     specified explicitly in the stimulus design with BoxCar basis functions
    
    
    properties
        bandwidth;
    end
    methods
        function obj = nonlin_GLM( prevJob )
            if nargin > 0, obj.prevJob = prevJob; end;
            
            obj.name = 'nonlinear GLM with AR(P) filtering';
            obj.basis('default') = nirs.design.basis.nonlinearHRF();
            obj.bandwidth = .3;
        end
        
        function S = runThis( obj, data )
            
             vec = @(x) x(:);
            
            for i = 1:length(data)
                % get data
                d=data(i).data;
                Fs=data(i).Fs;
                
                probe = data(i).probe;
                
                % make sure data is in order
                [probe.link, idx] = sortrows(probe.link, {'source', 'detector', 'type'});
                
                basis=obj.basis;
                
                cnt=0;
                param0=[];
                for idx=1:basis.count
                    base=basis.values{idx};
                    x=base.invtform_params;
                    lst{idx}=cnt+[1:length(x)];
                    cnt=cnt+length(x);
                    param0=[param0; x];
                end
                
                tHRF=[0:1/Fs:18];
                % Temporarirly change the basis to the FIR model to get the
                % design matrix
                obj.basis=Dictionary();
                FIRmodel=nirs.design.basis.FIR;
                FIRmodel.binwidth=1;
                FIRmodel.nbins=length(tHRF);
                FIRmodel.isIRF=true;
                
                obj.basis('default')=FIRmodel;
                % get experiment design
                [X, ~] = obj.createX( data(i) );
                % Put the proper basis back
                obj.basis=basis;
                
                
                C = obj.getTrendMatrix( data(i).time );
                
                % check model
                obj.checkRank( [X C] )
                obj.checkCondition( [X C] )
                
                names=unique(nirs.getStimNames(data(i)));
                
                H={};
                for j=1:length(names)
                    if(basis.iskey(names{j}))
                        basefcn=@(x)assignbasis(basis(names{j}),x);
                    else
                        basefcn=@(x)assignbasis(basis('default'),x);
                    end
                    H{end+1}=@(x)convert(basefcn(x(lst{j})),[1; zeros(length(tHRF)-1,1)],tHRF);
                end
                fcn=@(beta,x)myfcn(beta,x,H,cnt);
                
                x = zeros(length(tHRF)*length(names)*size(d,2),length(names)*size(d,2));
                for i=1:length(names)*size(d,2)
                    x([1:length(tHRF)]+(i-1)*length(tHRF),i)=1;
                end
                
                xo=zeros(cnt+size(x,2),1);
                for j=1:length(names)
                    xo(lst{j})=param0;
                end
               
                b=zeros(size(x,1),1);
                dd=d;
                
                options=statset('fitnlm');
                options.Robust='on';
                options.MaxIter=50;
                for iter=1:4
                    disp(['Iteration ' num2str(iter) ' of 4']); 
                    stats = nirs.math.ar_irls( dd, [X C], round(4*Fs) );
                    b2=stats.beta(1:size(X,2),:);
                    
                    if(obj.bandwidth*2/Fs<1)
                        [fa,fb]=butter(4,obj.bandwidth*2/Fs);
                        b2=filtfilt(fa,fb,b2);
                    end
                    
                    b=b+reshape(b2,[],1);
                    hrf=[];
                    for i=1:length(H)
                        hrf=[hrf; H{i}(xo)];
                    end
                    stats2 = nirs.math.ar_irls( d, [X*hrf C], round(4*Fs) );
                    xo(cnt+1:end)=reshape(stats2.beta(1:length(H),:),[],1);
                    
                    warning('off','stats:nlinfit:ModelConstantWRTParam');
                    warning('off','stats:nlinfit:IterationLimitExceeded');
                    %nlm = fitnlm(x,b,fcn,xo,'Options',options);
                    coef=nirs.math.NonLinearModel.nlinfit(x,b,fcn,xo,'Options',options);
                    
                    hrf=[];
                    for i=1:length(H)
                        %hrf=[hrf; H{i}(nlm.Coefficients.Estimate)];
                        hrf=[hrf; H{i}(coef)];
                    end
                    stats = nirs.math.ar_irls( d, [X*hrf C], round(4*Fs) );
                    dd=d-[X*hrf C]*stats.beta;
                end
                
                                
                % put stats
                ncond = length(names);
                nchan = size(data(i).probe.link, 1);
                
                link = repmat( probe.link, [ncond 1] );
                cond = repmat(names(:)', [nchan 1]);
                cond = cond(:);
                
                S(i) = nirs.core.ChannelStats();
                
                S(i).variables = [link table(cond)];
                S(i).beta = vec( stats.beta(1:ncond,:)' );
                
                covb = zeros( nchan*ncond );
                for j = 1:nchan
                    idx = (0:ncond-1)*nchan + j;
                    covb(idx, idx) = stats.covb(1:ncond, 1:ncond, j);
                end
                
                S(i).covb = covb;
                
                S(i).dfe  = stats.dfe(1);
                
                S(i).description = data(i).description;
                
                S(i).demographics   = data(i).demographics;
                S(i).probe          = data(i).probe;
                
                bb=Dictionary;
                for idx=1:basis.count
                    bb(basis.keys{idx})=assignbasis(basis('default'),coef(lst{idx}));    
                end
                for j=1:length(names)
                    if(basis.iskey(names{j}))
                        basefcn=@(x)assignbasis(basis(names{j}),x);
                    else
                        basefcn=@(x)assignbasis(basis('default'),x);
                    end
                    H{end+1}=@(x)convert(basefcn(x(lst{j})),[1; zeros(length(tHRF)-1,1)],tHRF);
                end
                
                 
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
                
                S(i).basis.base=bb;
                S(i).basis.Fs=Fs;
                S(i).basis.stim=stim;
                
                
                % print progress
                obj.printProgress( i, length(data) )
            end
            
        end
        
    end
    
end

function d = myfcn(beta,x,H,cnt)
    d=[];
    for i=1:length(H)
        d=[d; H{i}(beta)];
    end
    d=repmat(d,size(x,2)/length(H),1);
    d = d .*(x*beta(cnt+1:end));
    
end


function basis = assignbasis(basis,beta)
 
    basis=basis.fwdtform_params(beta);

end