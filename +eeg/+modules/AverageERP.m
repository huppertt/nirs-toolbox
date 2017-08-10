classdef AverageERP < eeg.modules.AbstractGLM
%% Average ERP - Performs first-level averaging of the ERP response

    properties
        time_window;   % time window to use for averaging
        prewhiten;
        ARorder;
        basis;
    end

    methods
        function obj = AverageERP( prevJob )
            if nargin > 0, obj.prevJob = prevJob; end
            
            obj.name = 'ERP Average';
            obj.time_window = [-0.050 0.300];
            obj.prewhiten=false;
            obj.ARorder = 10;
            obj.basis=nirs.design.basis.FIR;
        end
        
        function S = runThis( obj, data )
                      vec = @(x) x(:);
                      
            for i = 1:length(data)
                % get data
                d  = data(i).data;
                t  = data(i).time;
                Fs = data(i).Fs;
                
                probe = data(i).probe;
                
                basis=obj.basis;
                if(isa(basis,'nirs.design.basis.FIR'))
                    nWin = (obj.time_window(2)-obj.time_window(1))*data(i).Fs;
                    basis.binwidth=1;
                    basis.nbins=nWin;
                    basis.isIRF=false;
                    nWin=max(-obj.time_window(1),0)*data(i).Fs;
                else
                    nWin=[];
                end
                model=Dictionary;
                model('default')=basis;
                
                % get experiment design
                [X, names] = obj.createX( data(i),model);
                C = obj.getTrendMatrix( t );
                
                X(1:nWin,:)=[];
                C(1:nWin,:)=[];
                X=[X; zeros(nWin,size(X,2))];
                C=[C; zeros(nWin,size(C,2))];
               
                
                % check model
                obj.checkRank( [X C] )
                obj.checkCondition( [X C] )
                
                % run regression
                if(obj.prewhiten)
                    stats = nirs.math.ar_irls( d, [X C], obj.ARorder);
                else
                    stats.beta = [X C]\ d;
                    
                    iX=pinv([X C]'*[X C]);
                    for j = 1:size(d,2)
                        stats.covb(:,:,j) =  iX * var(d(:,j) - [X C]*stats.beta(:,j));
                    end
                    stats.dfe = size(d,1) - rank([X C]);
                end
                
                if(~isempty(nWin))
                    stats.beta(1:size(X,2),:)=stats.beta(1:size(X,2),:)-ones(size(X,2),1)*mean(stats.beta(1:nWin,:));
                end
                
                % put stats
                ncond = length(names);
                nchan = size(data(i).probe.link, 1);
                
                electrodes = repmat( probe.link, [ncond 1] );
                cond = repmat(names(:)', [nchan 1]);
                cond = cond(:);
                
                S(i) = eeg.core.ChannelStats();
                
                S(i).variables = [electrodes table(cond)];
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
                
                stim=Dictionary;
                for j=1:data(i).stimulus.count;
                    ss=data(i).stimulus.values{j};
                    if(isa(ss,'nirs.design.StimulusEvents'))
                        s=nirs.design.StimulusEvents;
                        s.name=ss.name;
                        s.dur=mean(ss.dur);
                        s.onset=-obj.time_window(1);
                        s.amp=1;
                        stim(data(i).stimulus.keys{j})=s;
                    end
                end
                
                S(i).basis.base=basis;
                S(i).basis.Fs=Fs;
                S(i).basis.stim=stim;
                
                % print progress
                obj.printProgress( i, length(data) )
            end

        end
        
    end
    
end