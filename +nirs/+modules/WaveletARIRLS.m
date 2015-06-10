classdef WaveletARIRLS < nirs.modules.AbstractGLM
   
    methods
        function obj = WaveletARIRLS( prevJob )
            if nargin > 0, obj.prevJob = prevJob; end
            
            obj.name = 'Solve GLM in Wavelet Domain';
            obj.basis('default') = nirs.design.basis.Canonical();
        end
        
        function S = runThis( obj, data )
            for i = 1:length(data)
                % get data
                d  = data(i).data;
                t  = data(i).time;
                Fs = data(i).Fs;
                
                % get experiment design
                [X, names] = obj.createX( data(i) );
                C = obj.getTrendMatrix( t );
                
                % check model
                obj.checkRank( [X C] )
                obj.checkCondition( [X C] )
                
                % run regression
                clear thisS;
                for j = 1:size(d,2)                    
                    y = d(:,j);
                    
                    lev = wmaxlev(length(y), 'sym8');
                    [ywav, l] = wavedec(y, lev, 'sym8');
                    
                    tmp = [X C];
                    for k = 1:size(tmp, 2)
                        [Xwav(:,k), l] = wavedec(tmp(:,k), lev, 'sym8');
                    end
                    L = cumsum([0; l(1:lev+1)]);
                    
                    b = Xwav \ ywav;
                    
                    for iter = 1:10
                        rwav = ywav - Xwav*b;
                        
                        % thresholding
                        yf = zeros(size(ywav));
                        Xf = zeros(size(Xwav));
                        for k = 1:lev+1
                            idx = L(k)+1:L(k+1);
                            
                            [a, ares] = ar_fit(rwav(idx), 2 );
                            
                            yf(idx) = filter( [1; -a(2:end)], 1, ywav(idx) ) / (mad(ares,0)/0.6745);
                            Xf(idx,:) = filter( [1; -a(2:end)], 1, Xwav(idx,:) ) / (mad(ares,0)/0.6745);
                        end
                        
                        [b, s] = robustfit(Xf, yf, [], [], 'off');
                    end

                    thisS.beta(:,j)     = b;
                    thisS.covb(:,:,j)   = s.covb;
                    thisS.dfe           = s.dfe;
                end
                
                % put stats
                ncond = length(names);
                S(i) = nirs.core.ChannelStats();
                
                S(i).beta = thisS.beta(1:ncond,:);
                S(i).covb = thisS.covb(1:ncond, 1:ncond, :);
                S(i).dfe  = thisS.dfe(1);
                S(i).description = data(i).description;
                
                S(i).names          = names';
                S(i).demographics   = data(i).demographics;
                S(i).probe          = data(i).probe;
                
                % print progress
                obj.printProgress( i, length(data) )
            end

        end
        
    end
    
end

