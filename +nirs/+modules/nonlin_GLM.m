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
    
    methods
        function obj = nonlin_GLM( prevJob )
            if nargin > 0, obj.prevJob = prevJob; end;
            
            obj.name = 'nonlinear GLM with AR(P) filtering';
            obj.basis('default') = nirs.design.basis.nonlinearHRF();
        end
        
        function S = runThis( obj, data )
            
            
            
            for i = 1:length(data)
                % get data
                d=data(i).data;
                Fs=data(i).Fs;
                
                probe = data(i).probe;
                
                % make sure data is in order
                [probe.link, idx] = sortrows(probe.link, {'source', 'detector', 'type'});
                
                basis=obj.basis;
                
                cnt=1;
                param0=[];
                for idx=1:basis.count
                    base=basis.values{idx};
                    x=base.invtform_params;
                    lst{idx}=cnt+[1:length(x)];
                    cnt=cnt+length(x);
                    param0=[param0; x];
                end
                
                
               
                % Temporarirly change the basis to the FIR model to get the
                % design matrix
                obj.basis=Dictionary();
                FIRmodel=nirs.design.basis.FIR;
                FIRmodel.binwidth=1;
                FIRmodel.nbins=1;
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
                
                
                HRFfcn={};
                XX={}; beta={}; Fcn={};
                for i=1:size(d,2)
                    H={};
                    for j=1:length(names)
                        if(basis.iskey(names{j}))
                            basefcn=@(x)assignbasis(basis(names{j}),x);
                        else
                            basefcn=@(x)assignbasis(basis('default'),x);
                        end
                        lst=cnt;
                        H{end+1}=@(x)x(lst)*basefcn(x);
                        cnt=cnt+1;
                    end
                    lst=cnt+[1:size(C,2)];
                    
                    H{end+1}=@(x)x(lst);
                    cnt=cnt+size(C,2);
                    
                    XX{end+1}=sparse([X C]);
                    beta{end+1}=d(:,i);
                    Fcn{end+1}=@(x)vertcat(cell2mat(cellfun(@(a){a(x)},H)'));
                end
                
                
                xo=zeros(cnt,1);
                xo(length(names)*size(d,2)+[1:length(param0)])=param0;
                
                hrf=@(x)vertcat(cell2mat(cellfun(@(a){a(x)},Fcn)'));
                
                options=statset('fitnlm');
                options.Robust='on';
                options.MaxIter=10;
                
                %             for iter=1:10
                %                 disp(iter)
                X=blkdiag(XX{:});
                b=vertcat(beta{:});
                
                model=@(beta,A)A*hrf(beta);
                [coef,r,J,Sigma,mse,errorModelInfo,robustw] = ...
                    nirs.math.NonLinearModel.nlinfit(X,b,model,xo,options);
                
                %                 nlm=fitnlm(X,b,@(beta,A)A*hrf(beta),xo,'Options',options);
                %
                %                 %xo=nlm.Coefficients.Estimate;
                %                 xo = coef;
                %
                %                 res = cellfun(@(b,X,F)(b-X*F(xo)),beta,XX,Fcn,'UniformOutput',false);
                %
                %                 pmax=Fs*4;
                %                 [~,f] = nirs.math.innovations(horzcat(res{:}),pmax);
                %
                %                 XX = cellfun(@(X,f)filter(f,1,full(X)),XX,f,'UniformOutput',false);
                %                 beta = cellfun(@(beta,f)filter(f,1,beta),beta,f,'UniformOutput',false);
                %
                %             end
                %
                
                % run regression
                stats = nirs.math.ar_irls( d, [X C], round(4*Fs) );
                
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
                
                % print progress
                obj.printProgress( i, length(data) )
            end
            
        end
        
    end
    
end


function J = jacobian(beta,obj,data,delta)

numparam=length(beta);
beta0=beta;

s2=getsigma2(obj,data);

for i=1:numparam
    beta=beta0;
    beta(i)=beta(i)+delta(i);
    obj.basis = assignbasis(obj.basis,beta);
    sigma2=getsigma2(obj,data);
    J(:,i)=(s2-sigma2)./delta(i);
end
end

function sigma2=getsigma2(obj,data)

d  = data.data;
time  = data.time;
Fs = data.Fs;

C = obj.getTrendMatrix( time );
[X, ~] = obj.createX( data(1) );
stats = nirs.math.ar_irls( d, [X C], round(4*Fs) );
sigma2=stats.sigma2;
end

function basis = assignbasis(basis,beta)
% 
% 
% cnt=0;
% for idx=1:basis.count
%     base=basis.values{idx};
%     x=base.invtform_params;
%     lst=cnt+[1:length(x)];
%     cnt=cnt+length(x);
%     
%     base=base.fwdtform_params(beta(lst));
%     basis(basis.keys{idx})=base;
% end
% 
% 
    basis=basis.fwdtform_params(beta);

end