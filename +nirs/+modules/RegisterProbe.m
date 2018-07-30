classdef RegisterProbe < nirs.modules.AbstractModule
%% RegisterProbe - Perform probe registration using supplied optode-reference and source-detector distances
%
% Usage:
% % Create probe registration job
% job = nirs.modules.RegisterProbe();
%
% % Supply distances between optodes and reference 10-20 points
% job.optode_reference_distance(end+1,:) = cell2table({'s1','Cz',1.2});
% job.optode_reference_distance(end+1,:) = cell2table({'s1','Fp1',2.3});
% job.optode_reference_distance(end+1,:) = cell2table({'d1','Fp1',2.5});
% job.optode_reference_distance(end+1,:) = cell2table({'d1','Fp2',5.2});
%
% % Supply distances between source-detector pairs
% job.source_detector_distance(end+1,:) = cell2table({'s1','d1',3.2});
% job.source_detector_distance(end+1,:) = cell2table({'s1','d2',3.5});
% job.source_detector_distance(end+1,:) = cell2table({'s2','d1',2.8});
% job.source_detector_distance(end+1,:) = cell2table({'s2','d2',3.7});
% % OR
% job.source_detector_distance = 3; % If the source-detector distance is fixed for all channels
%
% % Run the job
% hb_registered = job.run( hb_unregistered );
%
% Options: 
%     optode_reference_distance - [# x 3] table of distances between optodes and 10-20 references 
%                                 (optode, reference, distance)
%
%     source_detectior_distance - [#channels x 3] or scalar, distances between sources and detectors of each channel
%                                 (source, detector, distance)
%
%     units                     - 'cm' (default) or 'mm'

    properties
        optode_reference_distance = table({},{},[],'VariableNames',{'optode','reference','distance'});
        source_detector_distance = table({},{},[],'VariableNames',{'source','detector','distance'});
        units = 'cm';
    end
    
    methods
        function obj = RegisterProbe( prevJob )
         	if nargin > 0, obj.prevJob = prevJob; end
            obj.name = 'Register Probe';
        end
        
        function data = runThis( obj , data )
           
            % Check that all probes match the first one
            probe = data(1).probe;
            if numel(data)>1
                mismatch = false;
                for i = 2:numel(data)
                    if ~isequal( data(i).probe , probe )
                        mismatch = true;
                        break;
                    end
                end
                if mismatch, warning('Not all probes are the same. Using only the first for all files.'); end
            end
            
            % Perform registration
            probe = obj.registerProbe( probe );
            
            % Copy registered probe to all files
            for i = 1:length(data)
                data(i).probe = probe;
            end
            
        end
        
        function probe1020 = registerProbe( obj, probe )
            
            % Get Colin27 mesh and reference locations
            assert(isnumeric(probe.link.type),'Probe registration must be done before Beer-Lambert Law.');
            lambda = unique(probe.link.type,'stable');
            fwdBEM = nirs.registration.Colin27.BEM(lambda);
            fwdBEM.mesh(1).fiducials.Draw(:)=false;
            refs = fwdBEM.mesh(1).fiducials;

            % Get optode and constraint info
            optodes = probe.optodes;
            optodes = optodes( strcmp(optodes.Type,'Source') | strcmp(optodes.Type,'Detector'),:);
            num_optodes = height(optodes);
            
            % Fix optode names in Optode-Reference distance table
            num_ref_cons = height(obj.optode_reference_distance);
            for i = 1:num_ref_cons
                name = obj.optode_reference_distance.optode{i};
                if length(name)<6
                    num = str2double(name(2:end));
                    if strcmpi(name(1),'s')
                        type = 'Source';
                    elseif strcmpi(name(1),'d')
                        type = 'Detector';
                    else
                        error('Unknown type: %s',name);
                    end
                    obj.optode_reference_distance.optode{i} = sprintf('%s-%04i',type,num);
                end
            end
            
            % Populate Source-detector distances if distance is fixed
            if isscalar(obj.source_detector_distance)
                link = probe.link( probe.link.type==probe.link.type(1) , 1:2);
                dist = repmat(obj.source_detector_distance,height(link),1);
                obj.source_detector_distance = table(link.source,link.detector,dist,'VariableNames',{'source','detector','distance'});
            end
            
            % Fix source-detector names
            num_sd_cons = height(obj.source_detector_distance);
            [sources,detectors] = deal(cell(num_sd_cons,1));
            for i = 1:num_sd_cons
                source = obj.source_detector_distance.source(i);
                detector = obj.source_detector_distance.detector(i);
                if isnumeric(source)
                    source = sprintf('Source-%04i',source);
                elseif length(source{1})<6
                    num = str2double(source{1}(2:end));
                    source = sprintf('Source-%04i',num);
                end
                if isnumeric(detector)
                    detector = sprintf('Detector-%04i',detector);
                 elseif length(detector{1})<6
                    num = str2double(detector{1}(2:end));
                    detector = sprintf('Detector-%04i',num);
                end
                sources{i} = source;
                detectors{i} = detector;
            end
            obj.source_detector_distance.source = sources;
            obj.source_detector_distance.detector = detectors;
            
            % Fix units
            if strcmpi( obj.units , 'cm' )
                obj.optode_reference_distance.distance = 10 * obj.optode_reference_distance.distance;
                obj.source_detector_distance.distance = 10 * obj.source_detector_distance.distance;
            else
                assert(strcmpi( obj.units , 'mm' ),'Unrecognized units: %s',obj.units);
            end
            if all(strcmpi(refs.Units,'cm'))
                refs(:,2:4) = array2table( 10 * table2array(refs(:,2:4)));
            else
                assert(all(strcmpi(refs.Units,'mm')),'Inconsistent reference units');
            end
            
            % Create Optode-Reference constraint matrices
            % norm(optode_projmat * optode_coords - ref_locs) ~= optode_ref_dists
            optode_projmat = zeros(num_ref_cons,num_optodes); % Projection to replicate optode coordinates for each constraint
            ref_locs = zeros(num_ref_cons,3); % XYZ coordinats of reference points
            optode_ref_dists = zeros(num_ref_cons,1); % Distance between optode and reference
            for i = 1:num_ref_cons
                optode_ind = strcmpi( optodes.Name , obj.optode_reference_distance.optode{i} );
                reference_ind = strcmpi( refs.Name , obj.optode_reference_distance.reference{i} );
                optode_projmat(i,optode_ind) = 1;
                ref_locs(i,:) = table2array(refs(reference_ind,2:4));
                optode_ref_dists(i) = obj.optode_reference_distance.distance(i);
            end
            
            % Create Source-Detector distance constraint matrices
            % norm(source_projmat * optode_coords - detector_projmat * optode_coords) ~= sd_dists
            source_projmat = zeros(num_sd_cons,num_optodes);
            detector_projmat = zeros(num_sd_cons,num_optodes);
            sd_dists = nan(num_sd_cons,1);
            for i = 1:num_sd_cons
                source_ind = strcmpi( optodes.Name , obj.source_detector_distance.source{i} );
                detector_ind = strcmpi( optodes.Name , obj.source_detector_distance.detector{i} );
                source_projmat(i,source_ind) = 1;
                detector_projmat(i,detector_ind) = 1;
                sd_dists(i) = obj.source_detector_distance.distance(i);
            end
            
            % Construct matrix to project source-detector distances back to optode-space
            channel_invmat = pinv(source_projmat) + pinv(detector_projmat);
            
            % Convert mesh to logical volume, and get voxel coordinates
            % (to simplify point-to-mesh distance calculation)
            FV = struct;
            FV.faces = fwdBEM.mesh(1).faces;
            FV.vertices = fwdBEM.mesh(1).nodes;
            img = nirs.modules.RegisterProbe.polygon2voxel( FV , [181 217 181] , 'center' , false);
            [I,J,K] = ind2sub( [181 217 181] , find(img) );
            grid_coords = [I J K];
            grid_coords = bsxfun( @minus , grid_coords , [181 217 181]./2 );
            grid_coords = permute(grid_coords,[3 2 1]);

            % Find x that minimizes the sum of absolute differences between the
            % current optode-reference distances and those supplied
            f = @(x) sum(abs(sqrt(sum((optode_projmat * x - ref_locs).^2,2)) - optode_ref_dists));
            
            % Constrain x such that values must be <= 1 mm away from the
            % nearest scalp voxel
            C = @(x) min(sqrt(sum(bsxfun(@minus,x,grid_coords).^2,2)),[],3) - 1;
            
            % Constrain x such thats source-detector distances are equal to
            % those supplied
            Ceq = @(x) channel_invmat * ( sqrt(sum((source_projmat * x - detector_projmat * x).^2,2)) - sd_dists );
            
            % Wraper for nonlinear constraints
            nonlcon = @(x) deal( repmat(C(reshape(x,[],3)),3,1) , repmat(Ceq(reshape(x,[],3)),3,1) );

            % Set bounds and initial points
            lb = repmat(min(grid_coords,[],3),num_optodes,1);
            ub = repmat(max(grid_coords,[],3),num_optodes,1);
            x0 = pinv(optode_projmat) * ref_locs;
            
            % Run the optimizer
            fprintf('Running constrained nonlinear optimization\n');
            opt = optimset('fmincon');
            opt.MaxFunEvals = 1e6;
            opt.MaxIterations = 1000;
            opt.Display = 'iter';
            max_iter = 5; iter = 0;
            while 1
                
                [x,~,exitflag] = fmincon(f,x0,[],[],[],[],lb,ub,nonlcon,opt);
                
                if exitflag>0
                    break;
                end
                
                assert(iter<max_iter,'Could not converge');
                warning('Restarting optimization...');
            end
            
            optodes(:,2:4) = array2table(x);
            optodes.Units(:) = {'mm'};           

            % Assign outputs
            probe1020 = nirs.core.Probe1020;
            probe1020.optodes = probe.optodes;
            probe1020.optodes_registered = optodes;
            probe1020.link = probe.link;
            probe1020 = probe1020.register_mesh2probe(fwdBEM.mesh);
            probe1020.opticalproperties = nirs.media.tissues.brain(0.7, 50,lambda);
            probe1020.defaultdrawfcn = '3D mesh (frontal)';
            
        end
        
    end
    methods(Static)
        function Volume=polygon2voxel(FV,VolumeSize,mode,Yxz)
            % This function POLYGON2VOXEL will convert a Triangulated Mesh into a
            % Voxel Volume which will contain the discretized mesh. Discretization of a 
            % polygon is done by splitting/refining the face, until the longest edge
            % is smaller than 0.5 voxels. Followed by setting the voxel beneath the vertice 
            % coordinates of that small triangle to one.
            %
            % Volume=polygon2voxel(FV,VolumeSize,Mode,Yxz);
            %
            % Inputs,
            %   FV : A struct containing FV.faces with a facelist Nx3 and FV.vertices
            %        with a Nx3 vertice list. Such a structure is created by Matlab
            %        Patch function
            %   VolumeSize : The size of the output volume, example [100 100 100]
            %   Mode : (optional) if set to:
            %               'none', The vertices data is directly used as coordinates
            %                       in the voxel volume.
            %               'wrap', The vertices data is directly used as coordinates
            %                       in the voxel volume, coordinates outside are 
            %                       circular wrapped to the inside.
            %               'auto', The vertices data is translated and 
            %                       scaled with a scalar to fit inside the new volume.
            %               'center', coordinate 0,0,0 is set as the center of the volume
            %                       instead of the corner of the voxel volume.
            %               'clamp', The vertices data is directly used as coordinates
            %                       in the voxel volume, coordinates outside are 
            %                       clamped to the inside.
            %   (Optional)
            %   Yxz : If true (default) use Matlab convention 1th dimension Y, 
            %         2th dimension X, and last dimension Z. Otherwise 1th 
            %         dimension is X, 2th Y and last Z.
            %
            %                        
            % Outputs,
            %   Volume : The 3D logical volume, with all voxels part of the discretized
            %           mesh one, and all other voxels zero.
            %
            % Example,
            %   % Compile the c-coded function
            %   mex polygon2voxel_double.c -v
            %
            %   % Load a triangulated mesh of a sphere
            %   load sphere; 
            %
            %   % Show the mesh
            %   figure, patch(FV,'FaceColor',[1 0 0]); axis square;
            %
            %   % Convert the mesh to a voxelvolume
            %   Volume=polygon2voxel(FV,[50 50 50],'auto');
            %
            %   % Show x,y,z slices
            %   figure,
            %   subplot(1,3,1), imshow(squeeze(Volume(25,:,:)));
            %   subplot(1,3,2), imshow(squeeze(Volume(:,25,:)));
            %   subplot(1,3,3), imshow(squeeze(Volume(:,:,25)));
            %
            %   %  Show iso surface of result
            %   figure, patch(isosurface(Volume,0.1), 'Facecolor', [1 0 0]);
            %
            % Example2,
            %   % Compile the c-coded function
            %   mex polygon2voxel_double.c -v
            %
            %   % Make A Volume with a few blocks
            %   I = false(120,120,120);
            %   I(40:60,50:70,60:80)=1; I(60:90,45:75,60:90)=1;
            %   I(20:60,40:80,20:60)=1; I(60:110,35:85,10:60)=1;
            %
            %   % Convert the volume to a triangulated mesh
            %   FV = isosurface(I,0.8);
            %
            %   % Convert the triangulated mesh back to a surface in a volume
            %   J = polygon2voxel(FV,[120, 120, 120],'none'); 
            %   % Fill the volume
            %   J=imfill(J,'holes');
            % 
            %   % Differences between original and reconstructed
            %   VD = abs(J-I);
            %
            %   % Show the original Mesh and Mesh of new volume
            %   figure, 
            %   subplot(1,3,1),  title('original')
            %     patch(FV,'facecolor',[1 0 0],'edgecolor','none'), camlight;view(3);
            %   subplot(1,3,2), title('converted');
            %     patch(isosurface(J,0.8),'facecolor',[0 0 1],'edgecolor','none'), camlight;view(3);
            %   subplot(1,3,3), title('difference');
            %     patch(isosurface(VD,0.8),'facecolor',[0 0 1],'edgecolor','none'), camlight;view(3); 
            %
            % Function is written by D.Kroon University of Twente (May 2009)

            if(nargin<4), Yxz=true; end

            % Check VolumeSize size
            if(length(VolumeSize)==1)
                VolumeSize=[VolumeSize VolumeSize VolumeSize];
            end
            if(length(VolumeSize)~=3)
                error('polygon2voxel:inputs','VolumeSize must be a array of 3 elements ')
            end

            % Volume Size must always be an integer value
            VolumeSize=round(VolumeSize);

            sizev=size(FV.vertices);
            % Check size of vertice array
            if((sizev(2)~=3)||(length(sizev)~=2))
                error('polygon2voxel:inputs','The vertice list is not a m x 3 array')
            end

            sizef=size(FV.faces);
            % Check size of vertice array
            if((sizef(2)~=3)||(length(sizef)~=2))
                error('polygon2voxel:inputs','The vertice list is not a m x 3 array')
            end

            % Check if vertice indices exist
            if(max(FV.faces(:))>size(FV.vertices,1))
                error('polygon2voxel:inputs','The face list contains an undefined vertex index')
            end

            % Check if vertice indices exist
            if(min(FV.faces(:))<1)
                error('polygon2voxel:inputs','The face list contains an vertex index smaller then 1')
            end

            % Matlab dimension convention YXZ
            if(Yxz)
                FV.vertices=FV.vertices(:,[2 1 3]); 
            end

            switch(lower(mode(1:2)))
                case {'au'} % auto
                    % Make all vertices-coordinates positive
                    FV.vertices=FV.vertices-min(FV.vertices(:));
                    scaling=min((VolumeSize-1)./(max(FV.vertices(:))));
                    % Make the vertices-coordinates to range from 0 to 100
                    FV.vertices=FV.vertices*scaling+1;
                    Wrap=0;
                case {'ce'} % center
                    % Center the vertices
                    FV.vertices=FV.vertices+repmat((VolumeSize/2),size(FV.vertices,1),1);
                    Wrap=0;
                case {'wr'} %wrap
                    Wrap=1;
                case{'cl'} % clamp
                    Wrap=2;
                otherwise
                    Wrap=0;
            end

            % Separate the columns;
            FacesA=double(FV.faces(:,1));
            FacesB=double(FV.faces(:,2));
            FacesC=double(FV.faces(:,3));
            VerticesX=double(FV.vertices(:,1));
            VerticesY=double(FV.vertices(:,2));
            VerticesZ=double(FV.vertices(:,3));

            % Volume size to double
            VolumeSize=double(VolumeSize);

            % Call the mex function
            Volume=nirs.modules.RegisterProbe.polygon2voxel_double(FacesA,FacesB,FacesC,VerticesX,VerticesY,VerticesZ,VolumeSize,Wrap);
        end
        
        function Volume=polygon2voxel_double(FacesA,FacesB,FacesC,VerticesX,VerticesY,VerticesZ,VolumeSize,Wrap)
            Vertices=[VerticesX(:) VerticesY(:) VerticesZ(:)]-1;

            % List with all vertices coordinates of a face
            FaceVertices=[Vertices(FacesA,:) Vertices(FacesB,:) Vertices(FacesC,:)];
            Volume=false(VolumeSize);
            Volume=nirs.modules.RegisterProbe.DrawSplitFaces(FaceVertices,Volume,Wrap);
        end

        function Volume=DrawSplitFaces(FaceVertices,Volume,Wrap)
            VolumeSize=size(Volume);
            % Calculate squared edge distances
            dist1=(FaceVertices(:,1)-FaceVertices(:,4)).*(FaceVertices(:,1)-FaceVertices(:,4))+(FaceVertices(:,2)-FaceVertices(:,5)).*(FaceVertices(:,2)-FaceVertices(:,5))+(FaceVertices(:,3)-FaceVertices(:,6)).*(FaceVertices(:,3)-FaceVertices(:,6));
            dist2=(FaceVertices(:,7)-FaceVertices(:,4)).*(FaceVertices(:,7)-FaceVertices(:,4))+(FaceVertices(:,8)-FaceVertices(:,5)).*(FaceVertices(:,8)-FaceVertices(:,5))+(FaceVertices(:,9)-FaceVertices(:,6)).*(FaceVertices(:,9)-FaceVertices(:,6));
            dist3=(FaceVertices(:,1)-FaceVertices(:,7)).*(FaceVertices(:,1)-FaceVertices(:,7))+(FaceVertices(:,2)-FaceVertices(:,8)).*(FaceVertices(:,2)-FaceVertices(:,8))+(FaceVertices(:,3)-FaceVertices(:,9)).*(FaceVertices(:,3)-FaceVertices(:,9));

            % Calculate mFaceVertices(:,1) distance
            maxdist=max([dist1(:) dist2(:),dist3(:)],[],2);

            % Draw triangle if distance <=0.5 pixel
            check=maxdist>0.5;

            FVR=FaceVertices(~check,:);
            % Select Vertices which must be split
            FaceVertices=FaceVertices(check,:);
            if(~isempty(FaceVertices))
                dist1=dist1(check); 
                dist2=dist2(check); 
                dist3=dist3(check);

                DX=(FaceVertices(:,1)+FaceVertices(:,4))/2; DY=(FaceVertices(:,2)+FaceVertices(:,5))/2; DZ=(FaceVertices(:,3)+FaceVertices(:,6))/2;
                FA1=[DX,DY,DZ,FaceVertices(:,4),FaceVertices(:,5),FaceVertices(:,6),FaceVertices(:,7),FaceVertices(:,8),FaceVertices(:,9)];
                FB1=[FaceVertices(:,1),FaceVertices(:,2),FaceVertices(:,3),DX,DY,DZ,FaceVertices(:,7),FaceVertices(:,8),FaceVertices(:,9)];

                DX=(FaceVertices(:,1)+FaceVertices(:,7))/2; DY=(FaceVertices(:,2)+FaceVertices(:,8))/2; DZ=(FaceVertices(:,3)+FaceVertices(:,9))/2;
                FA2=[DX,DY,DZ,FaceVertices(:,4),FaceVertices(:,5),FaceVertices(:,6),FaceVertices(:,7),FaceVertices(:,8),FaceVertices(:,9)];
                FB2=[FaceVertices(:,1),FaceVertices(:,2),FaceVertices(:,3),FaceVertices(:,4),FaceVertices(:,5),FaceVertices(:,6),DX,DY,DZ];

                DX=(FaceVertices(:,7)+FaceVertices(:,4))/2; DY=(FaceVertices(:,8)+FaceVertices(:,5))/2; DZ=(FaceVertices(:,9)+FaceVertices(:,6))/2;
                FA3=[FaceVertices(:,1),FaceVertices(:,2),FaceVertices(:,3),DX,DY,DZ,FaceVertices(:,7),FaceVertices(:,8),FaceVertices(:,9)];
                FB3=[FaceVertices(:,1),FaceVertices(:,2),FaceVertices(:,3),FaceVertices(:,4),FaceVertices(:,5),FaceVertices(:,6),DX,DY,DZ];

                DX=(FaceVertices(:,1)+FaceVertices(:,7))/2; DY=(FaceVertices(:,2)+FaceVertices(:,8))/2; DZ=(FaceVertices(:,3)+FaceVertices(:,9))/2;
                FA4=[DX,DY,DZ,FaceVertices(:,4),FaceVertices(:,5),FaceVertices(:,6),FaceVertices(:,7),FaceVertices(:,8),FaceVertices(:,9)];
                FB4=[FaceVertices(:,1),FaceVertices(:,2),FaceVertices(:,3),FaceVertices(:,4),FaceVertices(:,5),FaceVertices(:,6),DX,DY,DZ];

                dist12=dist1>dist2;
                dist12n=~dist12;
                FA1(dist12n,:)=FA3(dist12n,:);
                FB1(dist12n,:)=FB3(dist12n,:);
                FA2(dist12n,:)=FA4(dist12n,:);
                FB2(dist12n,:)=FB4(dist12n,:);
                dist1(dist12n,:)=dist2(dist12n,:);

                dist13=dist1>dist3;
                dist13n=~dist13;
                FA1(dist13n,:)=FA2(dist13n,:);
                FB1(dist13n,:)=FB2(dist13n,:);

                FaceVertices=[FA1;FB1];

                % Split / Draw Vertices
                Volume=nirs.modules.RegisterProbe.DrawSplitFaces(FaceVertices,Volume,Wrap);
            end

            % Draw remaining faces
            FaceVertices=FVR;

            if(Wrap==0)
                % Calculate 1D volume indices
                indexA=nirs.modules.RegisterProbe.mindex3(floor(FaceVertices(:,1)+0.5),floor(FaceVertices(:,2)+0.5), floor(FaceVertices(:,3)+0.5), VolumeSize(1), VolumeSize(2), VolumeSize(3), Wrap);
                indexB=nirs.modules.RegisterProbe.mindex3(floor(FaceVertices(:,4)+0.5),floor(FaceVertices(:,5)+0.5), floor(FaceVertices(:,6)+0.5), VolumeSize(1), VolumeSize(2), VolumeSize(3), Wrap);
                indexC=nirs.modules.RegisterProbe.mindex3(floor(FaceVertices(:,7)+0.5),floor(FaceVertices(:,8)+0.5), floor(FaceVertices(:,9)+0.5), VolumeSize(1), VolumeSize(2), VolumeSize(3), Wrap);

                % Remove outside vertices
                checkA=(FaceVertices(:,1)<0)|(FaceVertices(:,2)<0)|(FaceVertices(:,3)<0)|(FaceVertices(:,1)>(VolumeSize(1)-1))|(FaceVertices(:,2)>(VolumeSize(2)-1))|(FaceVertices(:,3)>(VolumeSize(3)-1));
                checkB=(FaceVertices(:,4)<0)|(FaceVertices(:,5)<0)|(FaceVertices(:,6)<0)|(FaceVertices(:,4)>(VolumeSize(1)-1))|(FaceVertices(:,5)>(VolumeSize(2)-1))|(FaceVertices(:,6)>(VolumeSize(3)-1));
                checkC=(FaceVertices(:,7)<0)|(FaceVertices(:,8)<0)|(FaceVertices(:,9)<0)|(FaceVertices(:,7)>(VolumeSize(1)-1))|(FaceVertices(:,8)>(VolumeSize(2)-1))|(FaceVertices(:,9)>(VolumeSize(3)-1));
                indexA(checkA)=[];
                indexB(checkB)=[];
                indexC(checkC)=[];

                % Draw the vertices
                Volume(indexA)=true;
                Volume(indexB)=true;
                Volume(indexC)=true;
            else
                 Volume(nirs.modules.RegisterProbe.mindex3(floor(FaceVertices(:,1)+0.5),floor(FaceVertices(:,2)+0.5), floor(FaceVertices(:,3)+0.5), VolumeSize(1), VolumeSize(2), VolumeSize(3), Wrap))=1;
                 Volume(nirs.modules.RegisterProbe.mindex3(floor(FaceVertices(:,4)+0.5),floor(FaceVertices(:,5)+0.5), floor(FaceVertices(:,6)+0.5), VolumeSize(1), VolumeSize(2), VolumeSize(3), Wrap))=1;
                 Volume(nirs.modules.RegisterProbe.mindex3(floor(FaceVertices(:,7)+0.5),floor(FaceVertices(:,8)+0.5), floor(FaceVertices(:,9)+0.5), VolumeSize(1), VolumeSize(2), VolumeSize(3), Wrap))=1;
            end
        end

        function index=mindex3(x, y, z, sizx, sizy, sizz, Wrap) 
            if(Wrap==1)
                % Positive modules 
                x=mod(x,sizx);
                y=mod(y,sizy);
                z=mod(z,sizz);
            elseif(Wrap>1)
                % Clamp 
                x=max(x,0); x=min(x,sizx-1);
                y=max(y,0); y=min(y,sizy-1);
                z=max(z,0); z=min(z,sizz-1);
            end
            index=z*sizx*sizy+y*sizx+x;
            % matlab
            index=index+1;
        end

    end
end