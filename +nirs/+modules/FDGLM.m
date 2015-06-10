classdef FDGLM < nirs.modules.AbstractGLM
    properties
        fmin = 0.01;
        fmax = 0.5;
    end
    methods
        function obj = FDGLM( prevJob )
            if nargin > 0, obj.prevJob = prevJob; end
            
            obj.name = 'Solve GLM in Freq Domain';
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
%                 X = [X ones(size(X,1),1)];
%                 % check model
%                 obj.checkRank( [X C] )
%                 obj.checkCondition( [X C] )
                
                % run regression
                clear thisS;
                for j = 1:size(d,2)                    
                    y = d(:,j);
                    b = robustfit(dct(X), dct(y),[],[],'off');
                    
                    for iter = 1:5
                        r = y - X*b;
                        [a, ar] = ar_fit(r, round(4*Fs));
                        a1 = a(1);
                        a = a(2:end);
                        varE = mad(ar,0) / 0.6745; varE = varE*varE;

                        f = (0:length(t)-1)' / length(t) * Fs/2;
                        z = zeros(size(f));
                        for k = 1:length(a)
                            z = z + a(k)*exp(-2*pi*f*1i*k);
                        end
                        pwr = varE ./ ((1-z).*conj(1-z));

                        Xw = diag(sqrt(1./pwr))*dct(X);
                        yw = diag(sqrt(1./pwr))*dct(y);

                        lst = f > obj.fmin & f < obj.fmax;
                        Xw = Xw(lst,:);
                        yw = yw(lst,:);
                        
                        [b, s] = robustfit(Xw, yw, [], [], 'off');
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

