classdef ImageReconMFX < nirs.functional.AbstractModule
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
  
    properties
        formula = 'beta ~ cond*group + (1|subject)';
        fwdModel = nirs.HashTable(); % key is subject name or atlas
        transMtx;
        mesh;
    end
    
    methods

        function obj = ImageReconMFX( prevJob )
           obj.name = 'Image Recon w/ Random Effects';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function G = execute( obj, S )
            demo = nirs.functional.createDemographicsTable( S );
            
            % TODO
            [x, z, names] = nirs.functional. ...
                parseWilkinsonFormula( obj.formula, demo );
            
            x = sparse(x);
            z = sparse(z);
            
            
            %% SVDS; PER FWD MODEL
            
            %% MASK; ACROSS FWD MODELS
            
            %% ASSEMBLE X & Z
            X = []; Z = [];
            for i = 1:size(x,1)
                sname = demo.subject(1);
                
                if obj.fwdModel.iskey( sname )
                    thisX = obj.fwdModel( sname );
                elseif obj.fwdModel.iskey('default')
                    thisX = obj.fwdModel( 'default' );
                else
                    error(['No forward model for subject: ' sname '.'])
                end
                
                X = [X; kron(thisX(:,mask),x(i,:))];
                Z = [Z; kron(thisX(:,mask),z(i,:))];

                
            end
            
            %% FITTING
                
            %% PUT STATS
            
        end
        
        function options = getOptions( obj )
            options = [];
        end
           
        function obj = putOptions( obj, options )
        end
        
    end
    
end

