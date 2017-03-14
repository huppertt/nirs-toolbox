classdef NirfastBEM
    
    properties
        mesh; 
        probe;
        prop;
        preK={};
        Fm = 0;
    end
    
    methods
        %% Constructor
        function obj = NirfastBEM( mesh, prop, probe, Fm )
            if nargin > 0, obj.mesh = mesh; end
            if nargin > 1, obj.prop = prop; end
            if nargin > 2, obj.probe = probe; end
            if nargin > 3, obj.Fm = Fm; end
        end
        
        %% Methods
        meas = measurement( obj );
        mesh = getNirfastMeshes( obj )
        [J,meas] = jacobian( obj, type );
        [J,meas] = layeredJacobian( obj );
       
        function obj = set.probe(obj,probe)
            if(~isa(probe,'nirs.core.Probe1020'))
                warning('probe must be a 3D registered probe');
            end
            if(all(probe.optodes.Z==0))
                disp('warning: changing probe to 3D using "swap_reg" function');
                probe=probe.swap_reg;
            end
            obj.probe=probe;
        end
        
        
        function meshBEM=mesh2BEM(obj);
            meshBEM=nirs.core.Mesh;
            
            for idx=1:length(obj.mesh)    
                meshBEM.faces=[meshBEM.faces; obj.mesh(idx).faces+length(meshBEM.nodes)];
                meshBEM.nodes=[meshBEM.nodes; obj.mesh(idx).nodes];
            end
            
            meshBEM.regions=zeros(length(meshBEM.faces),2);
            cnt=0;
            for idx=1:length(obj.mesh)
               meshBEM.regions(cnt+[1:length(obj.mesh(idx).faces)],2)=idx;
               meshBEM.regions(cnt+[1:length(obj.mesh(idx).faces)],1)=idx-1;
               cnt=cnt+length(obj.mesh(idx).faces);
            end
         
            
        end
        
        function mesh = combinemesh(obj)
            mesh = nirs.core.Mesh;
        
           mesh=obj.mesh(1);
           
            for idx=2:length(obj.mesh)
                n=size(mesh.nodes,1);
               mesh.nodes=[mesh.nodes; obj.mesh(idx).nodes];
               mesh.faces=[mesh.faces; obj.mesh(idx).faces+n];
               mesh.elems=[mesh.elems; obj.mesh(idx).elems+n];
               mesh.regions=[mesh.regions; obj.mesh(idx).regions];
               try; mesh.fiducials=[mesh.fiducials; obj.mesh(idx).fiducials]; end;
            end
            
            
        end
        
        function obj = loadBEM_atlasviewer(obj,filename)
        
            atlas = load(filename);
            
            f=figure('visible','off');
            p=patch('vertices',atlas.headsurf.mesh.vertices,'faces',atlas.headsurf.mesh.faces);
            p=reducepatch(p,.2);
            atlas.headsurf.mesh.vertices=p.vertices;
            atlas.headsurf.mesh.faces=p.faces;
            
             p=patch('vertices',atlas.pialsurf.mesh.vertices,'faces',atlas.pialsurf.mesh.faces);
            p=reducepatch(p,.2);
            atlas.pialsurf.mesh.vertices=p.vertices;
            atlas.pialsurf.mesh.faces=p.faces;
            close(f);
            
            BEM(1)=nirs.core.Mesh();
            T = atlas.headsurf.T_2ras(1:3,1:3);
            BEM(1).nodes=atlas.headsurf.mesh.vertices;
            BEM(1).nodes=(BEM(1).nodes-ones(size(BEM(1).nodes,1),1)*atlas.headsurf.T_2ras(1:3,4)')*T;
            BEM(1).faces=atlas.headsurf.mesh.faces;
            
            Pos = atlas.refpts.pos;
            
            T = atlas.headsurf.T_2ras(1:3,1:3);
            Pos=(Pos-ones(size(Pos,1),1)*atlas.headsurf.T_2ras(1:3,4)')*T;

            marker{1}= atlas.refpts.labels';
           
            [TR, TT] = icp(BEM(1).nodes',Pos');
            Pos=(TR*Pos'+TT*ones(1,size(Pos,1)))';
%             k=dsearchn(BEM(1).nodes,Pos);
%             Pos=BEM(1).nodes(k,:); 
             
            fidtbl=table(marker{1},Pos(:,1),Pos(:,2),Pos(:,3),repmat({'10-20'},length(marker{1}),1),...
                repmat({'mm'},length(marker{1}),1),repmat(true,length(marker{1}),1),...
                'VariableNames',BEM(1).fiducials.Properties.VariableNames);
            
            if(height(BEM(1).fiducials)==0)
                BEM(1).fiducials=fidtbl;
            else
                BEM(1).fiducials=[BEM(1).fiducials; fidtbl];
            end
            
            BEM(2)=nirs.core.Mesh();
            T = atlas.headsurf.T_2ras(1:3,1:3);
            BEM(2).nodes=atlas.pialsurf.mesh.vertices;
            BEM(2).nodes=(BEM(2).nodes-ones(size(BEM(2).nodes,1),1)*atlas.headsurf.T_2ras(1:3,4)')*T;
            BEM(2).faces=atlas.pialsurf.mesh.faces;
            
            if(isfield(atlas,'probe'))
                atlas.probe.MeasList=atlas.probe.ml;
                
                optpos_reg=atlas.probe.optpos_reg;
                T = atlas.headsurf.T_2ras(1:3,1:3);
                optpos_reg=(optpos_reg-ones(size(optpos_reg,1),1)*atlas.headsurf.T_2ras(1:3,4)')*T;
                atlas.probe.DetPos=optpos_reg(1+size(atlas.probe.srcpos,1):...
                    size(atlas.probe.srcpos,1)+size(atlas.probe.detpos,1),:);
                atlas.probe.SrcPos=optpos_reg(1:size(atlas.probe.srcpos,1),:);
                
                if(~isfield(atlas.probe,'lambda'))
                    atlas.probe.lambda=[690 830];
                end
                
                atlas.probe.Lambda=atlas.probe.lambda;
                obj.probe = nirs.util.sd2probe(atlas.probe);
                
                BEM(1).fiducials.Draw(:)=false;
                BEM(1).fiducials=[BEM(1).fiducials;...
                    [obj.probe.optodes table(repmat(true,height(obj.probe.optodes),1),...
                    'VariableNames',{'Draw'})]];
                
                 lambda = unique(obj.probe.link.type);
            obj.prop{1}=nirs.media.tissues.skin(lambda);
            obj.prop{2}=nirs.media.tissues.brain(.70,50,lambda);

            end
           
            
            obj.mesh=BEM;
        end
        
        
        function obj = loadBEM_fif(obj,filename)
            
            % These are the levels of the icosohedral
            vertLength=[12 42 162 642 2562 10242];
            faceLength=[20 80 320 1280 5120 20480];
       
            surf=mne_read_bem_surfaces(filename);
            
            BEM(1)=nirs.core.Mesh();
            BEM(1).nodes=double(surf(1).rr)*1000;
            BEM(1).faces=double(surf(1).tris);
            
            
            fid=fopen(which('ext1020.sfp'),'r');
            marker=textscan(fid,'%s\t%d\t%d\t%d');
            fclose(fid);
            Pos=double([marker{2} marker{3} marker{4}]);
            
            Pos = icbm_spm2tal(Pos);
            
            [TR, TT] = icp(BEM(1).nodes',Pos');
            Pos=(TR*Pos'+TT*ones(1,size(Pos,1)))';
            k=dsearchn(BEM(1).nodes,Pos);
            Pos=BEM(1).nodes(k,:); 
             
            fidtbl=table(marker{1},Pos(:,1),Pos(:,2),Pos(:,3),repmat({'10-20'},length(marker{1}),1),...
                repmat({'mm'},length(marker{1}),1),repmat(true,length(marker{1}),1),...
                'VariableNames',BEM(1).fiducials.Properties.VariableNames);
            
            if(height(BEM(1).fiducials)==0)
                BEM(1).fiducials=fidtbl;
            else
                BEM(1).fiducials=[BEM(1).fiducials; fidtbl];
            end
            
            for i=2:length(surf)
                BEM(i)=nirs.core.Mesh();
                BEM(i).nodes=double(surf(i).rr)*1000;
                BEM(i).faces=double(surf(i).tris);
            end
            
            BEM(1).transparency=.1;
            BEM(2).transparency=.1;
            
            obj.mesh=BEM;
            
            
        end
        
        function obj = load_freesurfer(obj,subjectdir,subjid)
            file=dir(fullfile(subjectdir,subjid,'bem','*-bem.fif'));
            file=fullfile(subjectdir,subjid,'bem',file(1).name)
            obj=obj.loadBEM_fif(file);
            
            
        end
        
        
        function draw(obj,values,varargin)
        
            if(~exist('values'))
                for idx=1:length(obj.mesh)
                    obj.mesh(idx).draw();
                end
            else
                cnt=0;
                for idx=1:length(obj.mesh)
                    v=values(cnt+[1:length(obj.mesh(idx).nodes)]);
                    obj.mesh(idx).draw(v,varargin{:});
                    hold on;
                    cnt=cnt+length(obj.mesh(idx).nodes);
                end
            end
            
        end
        
    end
    

        
end

