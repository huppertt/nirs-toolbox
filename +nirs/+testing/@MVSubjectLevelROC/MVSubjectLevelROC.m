classdef MVSubjectLevelROC
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name    = 'Multivariate GLM ROC';
        
        truth
        tstat
        
        phbo
        phbr
        pboth
        
        pipeline
        
        niter           = 20;
        
        beta        	= 3;    % uM
        stimLength      = 5;    % sec
        stimSpacing     = 10;   % s
        
    end
    
    methods
        function obj = MVSubjectLevelROC( pipeline )
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
            %obj.tstat = zeros(obj.niter*nChan,1);
            obj.truth   = zeros(obj.niter*nChan/2,1);
            obj.phbo    = zeros(obj.niter*nChan/2,1);
            obj.phbr    = zeros(obj.niter*nChan/2,1);
            obj.pboth   = zeros(obj.niter*nChan/2,1);
            
            txt0 = '';
            for iter = 1:obj.niter
                % choose data file
                iData = randi(nData);
                d = hb(iData).data;
                t = hb(iData).time;
                
                % choose channels
                iChan = randperm(nChan/2, floor(nChan/4) + mod(iter,2));
                
                % random stim design
                [X, stim] = nirs.testing.randDesignMat( t, obj.stimLength, obj.stimSpacing );
                
                d(:,iChan) = bsxfun(@plus, d(:,iChan), X*obj.beta);           	% hbo
                d(:,iChan + nChan/2) = bsxfun(@plus,d(:,iChan), -X*obj.beta);   % hbr
                
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
                obj.phbo( (iter-1)*nChan/2+1:iter*nChan/2 ) = S.p(1,:);
                obj.phbr( (iter-1)*nChan/2+1:iter*nChan/2 ) = S.p(2,:);
                obj.pboth( (iter-1)*nChan/2+1:iter*nChan/2 ) = S.ftest([1 1]).p;
                obj.truth((iter-1)*nChan/2 + iChan) = 1;
                                
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

            [tp0, fp0, th0] = roc( obj.truth, obj.phbo );
            [tp1, fp1, th1] = roc( obj.truth, obj.phbr );
            [tp2, fp2, th2] = roc( obj.truth, obj.pboth );
            
            figure(h(1)), hold on
            plot( fp0, tp0, 'r' )
            plot( fp1, tp1, 'b' )
            plot( fp2, tp2, 'g' )
            plot(linspace(0,1), linspace(0,1),'r--')
            xlabel('False Postive Rate','FontSize',14)
            ylabel('True Positive Rate','FontSize',14)
            legend('hbo', 'hbr', 'joint', 'random','Location','SouthEast')
            
            figure(h(2)), hold on
            plot( th0, fp0, 'r' )
            plot( th1, fp1, 'b' )
            plot( th2, fp2, 'g' )
            plot(linspace(0,1), linspace(0,1),'r--')
            xlabel('p','FontSize',14)
            ylabel('False Positive Rate','FontSize',14)
            legend('hbo', 'hbr', 'joint', 'ideal','Location','SouthEast')
        end        
        
%         function options = getOptions( obj )
%             options = [];
%         end
%            
%         function obj = putOptions( obj, options )
%         end
    end
    
end

