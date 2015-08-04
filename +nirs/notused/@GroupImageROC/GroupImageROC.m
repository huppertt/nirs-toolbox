classdef TestImageMFX
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name        = 'ROC Test Image Recon w/ MFX';
        
        truth
        tstat
        p
        
        pipeline
        
        niter   	= 10;
        
        snr        	= 1;
        
    end
    
    methods
        function obj = TestImageMFX( pipeline )
           if nargin > 0
               obj.pipeline = pipeline;
           end
        end
        
        function obj = run( obj )
            % jacobian
            J = obj.pipeline.jacobian('default');
            
            % mesh
            mesh = obj.pipeline.mesh;
            nNode = size( mesh.nodes,1 );
            
            % sphere template
            load([fileparts(which('nirs.functional.modules.ImageReconMFX')) '/template.mat'])
            sphere.nodes = template.vertex{end}'*100;
            sphere.faces = template.face{end}';
            
            % mask
            mask = sum( abs(J.hbo),1 )';
            mask = mask > 0.1 * max(mask);
            lst = find(mask);
            
            % regression formula
            obj.pipeline.formula = 'beta ~ group + cont + (1|subject)';
            
            % preallocate
            obj.tstat = zeros(obj.niter*nNode*2,3);
            obj.truth = zeros(obj.niter*nNode*2,1);
            
            txt = '';
            for iter = 1:obj.niter
                
                % center of activation
                ipnt = lst( randi(length(lst)) );
                isph = max(1,mod(ipnt,length(sphere.nodes)));
                d = pwdist( sphere.nodes(isph,:), sphere.nodes );
                if ipnt > nNode
                    d = [1e6*ones(size(d)); d];
                else
                    d = [d; 1e6*ones(size(d))];
                end
                
                % activation
                b = zeros( size(mask) );
                b( d<15 ) = 1;
                b = [b; b]; % hbo + hbr
                
                % subject-level stats
                s = max( abs( [J.hbo J.hbr]*b ) ) / obj.snr;
                for i = 1:5
                    c = normrnd(1,0.25);
                    S(i).beta = ( [J.hbo J.hbr]*1*c*b )' + normrnd(0,s,1,size(J.hbo,1));
                    S(i).demographics = Dictionary({'subject','group','cont'},{num2str(i),'1',c});
                    S(i).stimulus = Dictionary({'test'},{[]});
                    S(i).covb = s^2*ones(1,1,size(J.hbo,1));
                end
                
                for i = 6:10
                    c = normrnd(1,0.25);
                    S(i).beta = ( [J.hbo J.hbr]*2*c*b )' + normrnd(0,s,1,size(J.hbo,1));
                    S(i).demographics = Dictionary({'subject','group','cont'},{num2str(i),'2',c});
                    S(i).stimulus = Dictionary({'test'},{[]});
                    S(i).covb = s^2*ones(1,1,size(J.hbo,1));
                end
                
                
                % image recon
                G = obj.pipeline.run( S );
                
                % put stats
                idx = (iter-1)*2*nNode + 1:iter*2*nNode;
                
                n = length(G.names);
                obj.truth(idx)     = b > 0;
                obj.tstat(idx,:)   = reshape( G.tstat, [length(G.tstat)/n n] );
                %obj.p(idx,:)       = tcdf(-obj.tstat(idx,:), G.dfe);
                obj.p(idx,:)      = reshape( G.p, [length(G.p)/n n] )
                
                txt0 = txt;
                txt = sprintf('Finished %6i of %6i.\n',iter, obj.niter);
%                 fprintf( repmat('\b',[1 length(txt0)]) );
                fprintf( txt );
            end
        end
        
        function draw( obj, h )
            if nargin < 2
                for i = 1:2*size(obj.p,2)
                    h(i) = figure;
                end
            end
            
            names = {'group1','difference','continuous'};
            
            % mask
            J = obj.pipeline.jacobian('default');
            mask = sum( abs(J.hbo),1 )';
            mask = mask > 0.01 * max(mask);
            
            lst = repmat(mask,[length(obj.truth)/length(mask) 1]);
%             lst = ones(size(mask)) > 0;
            for i = 1:size(obj.p,2)

                [tp, fp, th] = roc( obj.truth(lst), obj.p(lst,i) );

                figure(h(i)), hold on
                plot( fp, tp )
                plot(linspace(0,1), linspace(0,1),'r--')
                xlabel('False Postive Rate','FontSize',14)
                ylabel('True Positive Rate','FontSize',14)
                legend('Pipeline', 'Random','Location','SouthEast')
                title( names{i} )
            end
            
            for i = 1:size(obj.p,2)
                
                [tp, fp, th] = roc( obj.truth(lst), obj.p(lst,i) );
                
                figure(h(size(obj.p,2)+i)), hold on
                plot( th, fp )
                plot(linspace(0,1), linspace(0,1),'r--')
                xlabel('p','FontSize',14)
                ylabel('False Positive Rate','FontSize',14)
                legend('Pipeline', 'Ideal','Location','SouthEast')
                title( names{i} )
                
            end
        end        
        
%         function options = getOptions( obj )
%             options = [];
%         end
%            
%         function obj = putOptions( obj, options )
%         end
    end
    
end

