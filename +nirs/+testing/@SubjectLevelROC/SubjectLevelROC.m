classdef SubjectLevelROC
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name    = 'ROC Test Subject Level';
        
        truth
        tstat
        p
        
        pipeline
        
        niter           = 20;
        
        beta        	= 3;    % uM
        stimLength      = 5;    % sec
        stimSpacing     = 10;   % s
        
    end
    
    methods
        function obj = SubjectLevelROC( pipeline )
           if nargin > 0
               obj.pipeline = pipeline;
           end
        end
        
        function obj = run( obj, data )
            nData = length(data);
            nChan = size( data(1).data,2 );
            
            % covert to hb
            j = nirs.modules.OpticalDensity();
            j = nirs.modules.BeerLambertLaw( j );
            hb = j.run( data );
            
            % preallocate
            obj.tstat = zeros(obj.niter*nChan,1);
            obj.truth = zeros(obj.niter*nChan,1);
            
            txt0 = '';
            for iter = 1:obj.niter
                % choose data file
                iData = randi(nData);
                d = hb(iData).data;
                t = hb(iData).time;
                
                % choose channels
                iChan = randperm(nChan,nChan/2);
                
                % random stim design
                tmin = min( t ) + obj.stimLength;
                tmax = max( t ) - 2*obj.stimLength;
                
                nrnd = round( 2*(tmax-tmin)/obj.stimSpacing );
                dt = exprnd(obj.stimLength, [nrnd 1]);
                
                onset = tmin + cumsum([0; dt]);
                onset = onset( onset < tmax );
                
                dur = obj.stimLength * ones(size(onset));
                
                amp = ones(size(dur));
                
                stim = nirs.design.StimulusEvents();
                stim.amp = amp;
                stim.dur = dur;
                stim.onset = onset;
                
                stim = Dictionary({'roc'},{stim});
                
                % add to data
                basis = Dictionary({'default'}, {nirs.design.basis.Canonical()});
                X = nirs.design.createDesignMatrix( stim, t, basis );
                
                d(:,iChan) = bsxfun(@plus,d(:,iChan),X*obj.beta);
                
                % create data
                j = nirs.modules.InverseBeerLambert();
                tmp = hb(iData);
                tmp.data = d;
                tmp.stimulus = stim;
                
                tmp = j.run( tmp );
                
                % convert OD to raw
                tmp.data = exp( - tmp.data ) + 100;
                
                % run pipeline
                S = obj.pipeline.run( tmp );
                
                % put stats
                obj.tstat( (iter-1)*nChan+1:iter*nChan )    = S.tstat(1,:)';
                obj.p( (iter-1)*nChan+1:iter*nChan )        = tcdf(-S.tstat(1,:)', S.dfe);
                obj.truth( (iter-1)*nChan + iChan )         = 1;
                
                txt = sprintf('Finished %6i of %6i.\n',iter, obj.niter);
                fprintf( repmat('\b',[1 length(txt0)]) );
                fprintf( txt );
            end
        end
        
        function draw( obj, h )
            if nargin < 2
                h(1) = figure;
                h(2) = figure;
            end

            [tp, fp, th] = roc( obj.truth, obj.p );
            
            figure(h(1)), hold on
            plot( fp, tp )
            plot(linspace(0,1), linspace(0,1),'r--')
            xlabel('False Postive Rate','FontSize',14)
            ylabel('True Positive Rate','FontSize',14)
            legend('Pipeline', 'Random','Location','SouthEast')
            
            figure(h(2)), hold on
            plot( th, fp )
            plot(linspace(0,1), linspace(0,1),'r--')
            xlabel('p','FontSize',14)
            ylabel('False Positive Rate','FontSize',14)
            legend('Pipeline', 'Ideal','Location','SouthEast')
        end        
        
%         function options = getOptions( obj )
%             options = [];
%         end
%            
%         function obj = putOptions( obj, options )
%         end
    end
    
end

