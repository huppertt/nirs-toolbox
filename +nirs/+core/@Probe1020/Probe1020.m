classdef Probe1020 < nirs.core.Probe
    %% Probe1020.  This is a class derived from the core Probe class which adds
    % methods for registration, display, and co-registration based on the 10-20
    % landmarking

    properties
        optodes_registered;   % Registered version of optodes
        braindepth;  % depth of brain for modeling
        opticalproperties;  % optical properties for modeling

        defaultdrawfcn; % Drawing function to use

    end
    properties( Dependent = true )
        headcircum;  % circumference of head  % make dependent
        AP_arclength; % Anterior-Posterior arc length
        LR_arclength; % Left-right arclength
        srcPos3D      % nsrc x 3 array of source positions (mm)
        detPos3D      % ndet x 3 array of detector positions (mm)
    end

    properties (Access = private)
        zoom; % Flag to zoom in or show full 10-20 probe

        AP_distance;  % Distance from Oz to Fpz
        LR_distance;  % Distance from LPA to RPA
        IS_distance;  % Distance from center to Cz
        labels;  % labels of 10-20 points
        pts1020;  % 10-20 points
        mesh;

    end

    methods

        function obj = Probe1020(probe,headsize)

            if(nargin>1)
                tbl=nirs.util.register_headsize(headsize);
            else
                tbl=nirs.util.list_1020pts('?');
            end
            obj.pts1020=[tbl.X tbl.Y tbl.Z];
            obj.labels=tbl.Name;

            % Find the default arc lengths
            pt(1,:)=obj.pts1020(find(ismember(lower(obj.labels),'lpa')),:);
            pt(2,:)=obj.pts1020(find(ismember(lower(obj.labels),'rpa')),:);
            pt(3,:)=obj.pts1020(find(ismember(lower(obj.labels),'cz')),:);
            pt(4,:)=obj.pts1020(find(ismember(lower(obj.labels),'nas')),:);
            pt(5,:)=obj.pts1020(find(ismember(lower(obj.labels),'iz')),:);

            obj.AP_distance=norm(pt(4,:)-pt(5,:));
            obj.LR_distance=norm(pt(1,:)-pt(2,:));
            obj.IS_distance=norm(pt(3,:)-.5*(pt(1,:)-pt(2,:)));

            obj.braindepth = 10;  % Fix and cite

            if(nargin>0 && ~isempty(probe))
                obj=obj.registerprobe(probe);
                lambda=unique(obj.optodes.type);
            else
                lambda=[808];
            end

            obj=obj.set_mesh;

            obj.opticalproperties = nirs.media.tissues.brain(lambda,0.7, 50);

            obj.defaultdrawfcn='draw1020';
        end


        function srcPos = get.srcPos3D(obj)
            %% This function returns the src pos (in mm)

            optodes=obj.optodes_registered;
            lst=[];
            for i=1:height(optodes)
                if(isempty(optodes.Name{i}))
                    lst=[lst; i];
                end
            end
            optodes(lst,:)=[];

            tbl=sortrows(optodes,{'Type','Name'});
            lst=find(ismember(tbl.Type,'Source'));

            found=[];
            for i=1:length(lst)
                sIdx=str2num(tbl.Name{lst(i)}(strfind(tbl.Name{lst(i)},'-')+1:end));
                srcPos(sIdx,:)=[tbl.X(lst(i)) tbl.Y(lst(i)) tbl.Z(lst(i))];
                if(strcmp(tbl.Units(lst(i)),'cm'))
                    srcPos(sIdx,:)=srcPos(sIdx,:)*10;
                end
                found(sIdx)=1;
            end
            srcPos(~found,:)=NaN;
            %
            %             srcPos=[tbl.X(lst) tbl.Y(lst) tbl.Z(lst)];
            %
            %              %Convert to mm if needed
            %             lstCM=find(ismember(tbl.Units(lst),{'cm'}));
            %             srcPos(lstCM,:)=srcPos(lstCM,:)*10;
        end

        function detPos = get.detPos3D(obj)

            optodes=obj.optodes_registered;
            lst=[];
            for i=1:height(optodes)
                if(isempty(optodes.Name{i}))
                    lst=[lst; i];
                end
            end
            optodes(lst,:)=[];

            %% This function returns the det pos (in mm)
            tbl=sortrows(optodes,{'Type','Name'});
            lst=find(ismember(tbl.Type,'Detector'));

            found=[];
            for i=1:length(lst)
                dIdx=str2num(tbl.Name{lst(i)}(strfind(tbl.Name{lst(i)},'-')+1:end));
                detPos(dIdx,:)=[tbl.X(lst(i)) tbl.Y(lst(i)) tbl.Z(lst(i))];
                if(strcmp(tbl.Units(lst(i)),'cm'))
                    detPos(dIdx,:)=detPos(dIdx,:)*10;
                end
                found(dIdx)=1;
            end
            detPos(~found,:)=NaN;

            %             detPos=[tbl.X(lst) tbl.Y(lst) tbl.Z(lst)];
            %
            %             %Convert to mm if needed
            %             lstCM=find(ismember(tbl.Units(lst),{'cm'}));
            %             detPos(lstCM,:)=detPos(lstCM,:)*10;
        end

        function obj = apply_tform_mesh(obj,tform)
            m=obj.mesh;
            for i=1:length(m)
                n=m(i).nodes;
                n(:,4)=1;
                n=n*(tform);
                m(i).nodes=n(:,1:3);

                if(~isempty(m(i).fiducials))
                    n=[m(i).fiducials.X m(i).fiducials.Y m(i).fiducials.Z];
                    n(:,4)=1;
                    n=n*(tform);
                    m(i).fiducials.X=n(:,1);
                    m(i).fiducials.Y=n(:,2);
                    m(i).fiducials.Z=n(:,3);
                end

            end
            n=obj.pts1020;
            n(:,4)=1;
            n=n*(tform);
            obj.pts1020=n(:,1:3);
            obj.mesh=m;

            % Find the default arc lengths
            pt(1,:)=obj.pts1020(find(ismember(lower(obj.labels),'lpa')),:);
            pt(2,:)=obj.pts1020(find(ismember(lower(obj.labels),'rpa')),:);
            pt(3,:)=obj.pts1020(find(ismember(lower(obj.labels),'cz')),:);
            pt(4,:)=obj.pts1020(find(ismember(lower(obj.labels),'nas')),:);
            pt(5,:)=obj.pts1020(find(ismember(lower(obj.labels),'iz')),:);

            obj.AP_distance=norm(pt(4,:)-pt(5,:));
            obj.LR_distance=norm(pt(1,:)-pt(2,:));
            obj.IS_distance=norm(pt(3,:)-.5*(pt(1,:)-pt(2,:)));
        end


        function obj = register_mesh2probe(obj,mesh,noreg)
            % This reshapes the mesh to fit the current probe

            if(nargin<3)
                noreg=false;
            end

            tbl=table(obj.labels,obj.pts1020(:,1),obj.pts1020(:,2),obj.pts1020(:,3),...
                'VariableNames',{'Name','X','Y','Z'});

            if(~noreg)
                if(height(mesh(1).fiducials)>3)
                    T = nirs.registration.cp2tform(mesh(1).fiducials,tbl);
                else
                    T=eye(4);
                end
                for idx=1:length(mesh)
                    n=mesh(idx).nodes;
                    n(:,4)=1;
                    n=n*T;
                    mesh(idx).nodes=n(:,1:3);

                    if(~isempty(mesh(idx).fiducials))
                        p=[mesh(idx).fiducials.X mesh(idx).fiducials.Y mesh(idx).fiducials.Z];
                        p(:,4)=1;
                        p=p*T;
                        mesh(idx).fiducials.X=p(:,1);
                        mesh(idx).fiducials.Y=p(:,2);
                        mesh(idx).fiducials.Z=p(:,3);


                    end

                end
            end
            obj.mesh=mesh;


        end


        function obj=SetFiducialsVisibility(obj,flag)
            if(nargin<2)
                flag=false;
            end
            obj.mesh(1).fiducials.Draw(:)=flag;
        end


        function obj=RealignMesh(obj)

            pt(1,:)=obj.pts1020(find(ismember(lower(obj.labels),'lpa')),:);
            pt(2,:)=obj.pts1020(find(ismember(lower(obj.labels),'rpa')),:);
            pt(3,:)=obj.pts1020(find(ismember(lower(obj.labels),'cz')),:);
            pt(4,:)=obj.pts1020(find(ismember(lower(obj.labels),'nas')),:);
            pt(5,:)=obj.pts1020(find(ismember(lower(obj.labels),'iz')),:);

            a=obj.LR_distance/2;
            b=obj.IS_distance;
            c=obj.AP_distance/2;

            d=mean(pt([1 2 4 5],3));

            R = [-a 0 (pt(2,3)+pt(1,3))/2
                a 0 (pt(2,3)+pt(1,3))/2
                0 0 pt(3,3)
                0 c pt(4,3)
                0 -c pt(5,3)];
            R(:,4)=1;
            pt(:,4)=1;
            RR=pt\R;

            for i=1:length(obj.mesh)
                n=obj.mesh(i).nodes;
                n(:,4)=1;
                n=n*RR;
                obj.mesh(i).nodes=n(:,1:3);

                if(~isempty(obj.mesh(i).fiducials))
                    n=[obj.mesh(i).fiducials.X obj.mesh(i).fiducials.Y obj.mesh(i).fiducials.Z];
                    n(:,4)=1;
                    n=n*RR;
                    obj.mesh(i).fiducials.X=n(:,1);
                    obj.mesh(i).fiducials.Y=n(:,2);
                    obj.mesh(i).fiducials.Z=n(:,3);

                end
            end


            n=[obj.optodes_registered.X obj.optodes_registered.Y obj.optodes_registered.Z];
            n(:,4)=1;
            n=n*RR;
            obj.optodes_registered.X=n(:,1);
            obj.optodes_registered.Y=n(:,2);
            obj.optodes_registered.Z=n(:,3);


            n=obj.pts1020;
            n(:,4)=1;
            n=n*RR;
            obj.pts1020=n(:,1:3);

        end




        function headsize=get_headsize(obj)
            headsize=Dictionary();
            headsize('lpa-cz-rpa')=obj.LR_arclength;
            headsize('Iz-cz-nas')=obj.AP_arclength;
            headsize('circumference')=obj.headcircum;

        end
        function varargout=draw(obj,varargin)

            if(~isempty(obj.link))
                if(ismember('hyperscan',obj.link.Properties.VariableNames))
                    p=obj;
                    p.link=p.link(ismember(p.link.hyperscan,'A'),:);

                    S={}; D={};
                    for i=1:height(p.link)
                        if iscell(p.link.source(i))
                            source = p.link.source{i};
                            detector = p.link.detector{i};
                        else
                            source = p.link.source(i);
                            detector = p.link.detector(i);
                        end
                        for j=1:length(source)
                            s=['000' num2str(source(j))];
                            S{end+1}=['Source-' s(end-3:end)];
                            d=['000' num2str(detector(j))];
                            D{end+1}=['Detector-' d(end-3:end)];
                        end
                    end

                    p.optodes_registered=p.optodes_registered(ismember(p.optodes.Name,{S{:} D{:}}),:);
                    p.optodes=p.optodes(ismember(p.optodes.Name,{S{:} D{:}}),:);
                    obj = p;
                end
            end

            if(~isempty(strfind(obj.defaultdrawfcn,'zoom')))
                obj.zoom=true;
            else
                obj.zoom=false;
            end

            if(~isempty(strfind(obj.defaultdrawfcn,'label')))
                addlabels=true;
            else
                addlabels=false;
            end

            if ~isempty(strfind(lower(obj.defaultdrawfcn),'ball')) && length(varargin)>1
                obj.defaultdrawfcn = '3D Mesh';
            end

            if(~isempty(strfind(obj.defaultdrawfcn,'10-20 map')));
                l=draw1020interp(obj,varargin{:});
            elseif( ~isempty(strfind(lower(obj.defaultdrawfcn),'ball')) & ~isempty(strfind(lower(obj.defaultdrawfcn),'3d')))
                if(strfind(obj.defaultdrawfcn,'('))
                    v=obj.defaultdrawfcn(strfind(obj.defaultdrawfcn,'(')+1:end);
                    v=v(1:strfind(v,')')-1);
                else
                    v='frontal';
                end
                varg=cell(2,1);
                varg{2}=addlabels;
                for ii=1:length(varargin)
                    varg{ii}=varargin{ii};
                end
                varargin=varg;

                l=draw3d_ballandstick(obj,varargin{:});
                axis_handle=[];
                for i=1:length(varargin)
                    if(isa(varargin{i},'matlab.ui.control.UIAxes') | isa(varargin{i},'matlab.graphics.axis.Axes'))
                        axis_handle=varargin{i};
                    end
                end
                if(isempty(axis_handle))
                    axis_handle=gca;
                end


                mesh=obj.getmesh;
                h=mesh.draw([],[],[],[],axis_handle);
                axis(axis_handle,'tight');
                nirs.util.rotateview(get(h(1),'Parent'),v);

            elseif(~isempty(strfind(obj.defaultdrawfcn,'3D')));

                if(strfind(obj.defaultdrawfcn,'('))
                    v=obj.defaultdrawfcn(strfind(obj.defaultdrawfcn,'(')+1:end);
                    v=v(1:strfind(v,')')-1);
                else
                    v='frontal';
                end

                varg=cell(4,1);
                varg{4}=addlabels;
                for ii=1:length(varargin)
                    varg{ii}=varargin{ii};
                end
                varargin=varg;
                
                l=draw3d(obj,varargin{:});
                if(~isempty(strfind(obj.defaultdrawfcn,'mesh')))
                    axis_handle=[];
                    for i=1:length(varargin)
                        if(isa(varargin{i},'matlab.ui.control.UIAxes') | isa(varargin{i},'matlab.graphics.axis.Axes'))
                            axis_handle=varargin{i};
                        end
                    end
                    if(isempty(axis_handle))
                        axis_handle=gca;
                    end


                    mesh=obj.getmesh;
                    h=mesh.draw([],[],[],[],axis_handle);
                    axis(axis_handle,'tight');
                    nirs.util.rotateview(get(h(1),'Parent'),v);
                end
            elseif(~isempty(strfind(obj.defaultdrawfcn,'10-20')));
                varg=cell(4,1);
                varg{4}=addlabels;
                for ii=1:length(varargin)
                    varg{ii}=varargin{ii};
                end
                varargin=varg;
                l=draw1020(obj,varargin{:});
            else
                l=draw@nirs.core.Probe(obj,varargin{:});
            end

            if(nargout>0)
                varargout{1}=l;
            end
        end

        function str = get.defaultdrawfcn(obj)
            str = obj.defaultdrawfcn;
        end
        function obj = set.defaultdrawfcn(obj,str)

            if(nargin==1)
                return
            end

            allowed={'10-20','10-20 mercator projection map';...
                '10-20 zoom', '10-20 mercator with restricted view';...
                '10-20 map', '10-20 mercator with underlain image';...
                '10-20 map zoom', '10-20 mercator with underlain image';...
                '10-20 label','10-20 mercator projection map';...
                '10-20 map zoom label', '10-20 mercator with underlain image';...
                '3D', '3D line drawing ';...
                '3D mesh', '3D line drawing overlain on mesh';...
                '3D mesh (frontal)', '3D line drawing overlain on mesh';...
                '3D mesh (left)', '3D line drawing overlain on mesh';...
                '3D mesh (right)', '3D line drawing overlain on mesh';...
                '3D mesh (superior)', '3D line drawing overlain on mesh';...
                '3D mesh (top)', '3D line drawing overlain on mesh';...
                '3D mesh (posterior)', '3D line drawing overlain on mesh';...
                '3D mesh (occipital)', '3D line drawing overlain on mesh';...
                '3D mesh (back)', '3D line drawing overlain on mesh';...
                '3D mesh (inferior)', '3D line drawing overlain on mesh';...
                '3D mesh (bottom)', '3D line drawing overlain on mesh';...
                '3D ball', '3D ball and stick drawing overlain on mesh';...
                '3D label', '3D line drawing ';...
                '3D label mesh', '3D line drawing overlain on mesh';...
                '3D label mesh (frontal)', '3D line drawing overlain on mesh';...
                '3D labelmesh (left)', '3D line drawing overlain on mesh';...
                '3D label mesh (right)', '3D line drawing overlain on mesh';...
                '3D label mesh (superior)', '3D line drawing overlain on mesh';...
                '3D label mesh (top)', '3D line drawing overlain on mesh';...
                '3D label mesh (posterior)', '3D line drawing overlain on mesh';...
                '3D label mesh (occipital)', '3D line drawing overlain on mesh';...
                '3D label mesh (back)', '3D line drawing overlain on mesh';...
                '3D label mesh (inferior)', '3D line drawing overlain on mesh';...
                '3D label mesh (bottom)', '3D line drawing overlain on mesh';...
                '3D label ball', '3D ball and stick drawing overlain on mesh';...
                '2D', '2D probe layout'};

            if(~isempty(str))
               for ii=1:length(allowed)
                   allowedABC{ii}=strip(sort(allowed{ii,1}));
               end
               strABC=strip(sort(str));
              
               idx=find(ismember(lower({allowedABC{:}}),lower(strABC)));
            else
                idx=[];
            end

            if(~isempty(idx))
                obj.defaultdrawfcn=allowed{idx,1};
            elseif(isempty(str) || strcmp(str,'?'))
                disp('Here are the options for drawing configurations');
                disp(allowed);
            end
        end


        function obj=swap_reg(obj)
            op=obj.optodes;
            obj.optodes=obj.optodes_registered;
            obj.optodes_registered=op;
        end
        function mesh = getmesh(obj);
            % create the spherical mesh
            mesh=obj.mesh;
        end

        function obj=set_mesh(obj,mesh,tbl1020)

            if(nargin>1)
                obj.mesh=mesh;

                if(nargin>2)
                    obj.pts1020=[tbl1020.X tbl1020.Y tbl1020.Z];
                    obj.labels=tbl1020.Name;

                    % Find the default arc lengths
                    pt(1,:)=obj.pts1020(find(ismember(lower(obj.labels),'lpa')),:);
                    pt(2,:)=obj.pts1020(find(ismember(lower(obj.labels),'rpa')),:);
                    pt(3,:)=obj.pts1020(find(ismember(lower(obj.labels),'cz')),:);
                    pt(4,:)=obj.pts1020(find(ismember(lower(obj.labels),'nas')),:);
                    pt(5,:)=obj.pts1020(find(ismember(lower(obj.labels),'iz')),:);

                    obj.AP_distance=norm(pt(4,:)-pt(5,:));
                    obj.LR_distance=norm(pt(1,:)-pt(2,:));
                    obj.IS_distance=norm(pt(3,:)-.5*(pt(1,:)-pt(2,:)));
                end
            else

                mesh(1) = nirs.util.spheremesh;
                mesh(2) = nirs.util.spheremesh;
                mesh(3) = nirs.util.spheremesh;

                pt(1,:)=obj.pts1020(find(ismember(obj.labels,'nas')),:);
                pt(2,:)=obj.pts1020(find(ismember(obj.labels,'Iz')),:);
                com=.5*sum(pt,1);

                t1020=nirs.util.list_1020pts;
                t1020b=nirs.util.list_1020pts('?',obj.get_headsize);


                fract=.2;
                mesh(1).nodes(:,1)=mesh(1).nodes(:,1)*obj.LR_distance/2+com(1);
                mesh(1).nodes(:,2)=mesh(1).nodes(:,2)*obj.AP_distance/2+com(2);
                mesh(1).nodes(:,3)=mesh(1).nodes(:,3)*obj.IS_distance+com(3);

                mesh(2).nodes(:,1)=mesh(2).nodes(:,1)*(obj.LR_distance/2-obj.braindepth*fract)+com(1);
                mesh(2).nodes(:,2)=mesh(2).nodes(:,2)*(obj.AP_distance/2-obj.braindepth*fract)+com(2);
                mesh(2).nodes(:,3)=mesh(2).nodes(:,3)*(obj.IS_distance-obj.braindepth*fract)+com(3);

                mesh(3).nodes(:,1)=mesh(3).nodes(:,1)*(obj.LR_distance/2-obj.braindepth)+com(1);
                mesh(3).nodes(:,2)=mesh(3).nodes(:,2)*(obj.AP_distance/2-obj.braindepth)+com(2);
                mesh(3).nodes(:,3)=mesh(3).nodes(:,3)*(obj.IS_distance-obj.braindepth)+com(3);

                fidtbl=table(obj.labels,obj.pts1020(:,1),obj.pts1020(:,2),obj.pts1020(:,3),...
                    repmat({'10-20'},length(obj.labels),1),...
                    repmat({'mm'},length(obj.labels),1),...
                    repmat(true,length(obj.labels),1),...
                    'VariableNames',{'Name','X','Y','Z','Units','Type','Draw'});

                %             [TR, TT] = icp(mesh(1).nodes',obj.pts1020');
                %             mesh(1).nodes=(TR'*mesh(1).nodes'-TT*ones(1,size(mesh(1).nodes,1)))';
                %             mesh(2).nodes=(TR'*mesh(2).nodes'-TT*ones(1,size(mesh(2).nodes,1)))';
                %             mesh(3).nodes=(TR'*mesh(3).nodes'-TT*ones(1,size(mesh(3).nodes,1)))';

                for iter=1:20
                    [k,d]=dsearchn(mesh(1).nodes,obj.pts1020);
                    a=mesh(1).nodes(k,:);
                    b=obj.pts1020;

                    T=a\b;
                    for i=1:length(mesh)
                        n=mesh(i).nodes;
                        n=n*T;
                        mesh(i).nodes=n(:,1:3);
                    end
                end

                mesh(1).fiducials=fidtbl;
                mesh(1).transparency=.2;
                mesh(2).transparency=.1;
                obj.mesh=mesh;
            end

        end

        function varargout=draw3d(obj,colors, lineStyles, axis_handle,addlabels)

            if(~isempty(obj.link))
                if isnumeric(obj.link.type)
                    link = obj.link(obj.link.type==obj.link.type(1),1:2);
                else
                    link = obj.link(strcmp(obj.link.type,obj.link.type(1)),1:2);
                end
            else
                link=table;
            end

            n = height(link);

            if nargin < 2 || isempty(colors)
                colors = repmat([0.3 0.5 1], [n 1]);
            elseif size(colors,1) == 1
                colors = repmat(colors, [n 1]);
            end

            if nargin < 3 || isempty(lineStyles)
                lineStyles = repmat({'LineStyle', '-', 'LineWidth', 6}, [n 1]);
            elseif size(lineStyles, 1) == 1
                lineStyles = repmat({'LineStyle', '-', 'LineWidth', 6}, [n 1]);
            end

            if nargin < 4 || isempty(axis_handle)
                axis_handle = axes();
            end

            if(nargin<5)
                addlabels=false;
            end

            % Points from the probe
            Pos(:,1)=obj.optodes_registered.X;
            Pos(:,2)=obj.optodes_registered.Y;
            Pos(:,3)=obj.optodes_registered.Z;
            DetPos3D=obj.detPos3D;
            SrcPos3D=obj.srcPos3D;


            hold(axis_handle,'on');
            lstS=find(ismember(obj.optodes_registered.Type,'Source'));
            scatter3(axis_handle,Pos(lstS,1),Pos(lstS,2),Pos(lstS,3),'filled','MarkerFaceColor','r')
            lstD=find(ismember(obj.optodes_registered.Type,'Detector'));
            scatter3(axis_handle,Pos(lstD,1),Pos(lstD,2),Pos(lstD,3),'filled','MarkerFaceColor','b')


            if(addlabels)
                for i=1:length(lstS)
                    text(Pos(lstS(i),1),Pos(lstS(i),2),Pos(lstS(i),3),['S' num2str(i)],'HorizontalAlignment','center','parent',axis_handle, 'FontSize', 12);
                end
                for i=1:length(lstD)

                    text(Pos(lstD(i),1),Pos(lstD(i),2),Pos(lstD(i),3),['D' num2str(i)],'HorizontalAlignment','center','parent',axis_handle, 'FontSize', 12);
                end
            end




            h=[];
            for i=1:height(link)
                if iscell(link.source(i))
                    source = link.source{i};
                    detector = link.detector{i};
                else
                    source = link.source(i);
                    detector = link.detector(i);
                end
                for j=1:length(source)
                    s = source(j);
                    d = detector(j);
                    h(i)=line(axis_handle,[SrcPos3D(s,1) DetPos3D(d,1)],[SrcPos3D(s,2) DetPos3D(d,2)],[SrcPos3D(s,3) DetPos3D(d,3)],'Color', colors(i, :), lineStyles{i, :});
                    set(h(i),'UserData',[s d]);
                end
            end

            axis(axis_handle,'equal');

            if(nargout>0)
                varargout{1}=h;
            end

        end

        function varargout=draw3d_ballandstick(obj,axis_handle,addlabels)

            if(~isempty(obj.link))
                if isnumeric(obj.link.type)
                    link = obj.link(obj.link.type==obj.link.type(1),1:2);
                else
                    link = obj.link(strcmp(obj.link.type,obj.link.type(1)),1:2);
                end
            else
                link=table;
            end

            n = height(link);

            if nargin < 2 || isempty(axis_handle)
                axis_handle = axes();
            end
            if(nargin<3)
                addlabels=false;
            end

            % Points from the probe
            Pos(:,1)=obj.optodes_registered.X;
            Pos(:,2)=obj.optodes_registered.Y;
            Pos(:,3)=obj.optodes_registered.Z;
            DetPos3D=obj.detPos3D;
            SrcPos3D=obj.srcPos3D;

            hold(axis_handle,'on');

            % Create shapes for 3D probe rendering
            [x_sph,y_sph,z_sph] = sphere(50);
            x_sph = 4*x_sph;
            y_sph = 4*y_sph;
            z_sph = 4*z_sph;

            % Draw optodes
            lstS=find(ismember(obj.optodes_registered.Type,'Source'));
            lstD=find(ismember(obj.optodes_registered.Type,'Detector'));
            src_coord = Pos(lstS,:);
            det_coord = Pos(lstD,:);
            for i=1:size(src_coord,1)
                surf(axis_handle, x_sph+src_coord(i,1), y_sph+src_coord(i,2), z_sph+src_coord(i,3) ,[],'FaceColor',[1 0 0],'EdgeAlpha',0);
            end
            for i=1:size(det_coord,1)
                surf(axis_handle, x_sph+det_coord(i,1), y_sph+det_coord(i,2), z_sph+det_coord(i,3) ,[],'FaceColor',[0 0 1],'EdgeAlpha',0);
            end
            
            if(addlabels)
                Pos2=Pos+6*(Pos./(sqrt(Pos(:,1).^2+Pos(:,2).^2+Pos(:,3).^2)*ones(1,3)));
                for i=1:length(lstS)
                    text(Pos2(lstS(i),1),Pos2(lstS(i),2),Pos2(lstS(i),3),['S' num2str(i)],'HorizontalAlignment','center','parent',axis_handle, 'FontSize', 16);
                end
                for i=1:length(lstD)

                    text(Pos2(lstD(i),1),Pos2(lstD(i),2),Pos2(lstD(i),3),['D' num2str(i)],'HorizontalAlignment','center','parent',axis_handle, 'FontSize', 16);
                end
            end

            % Draw channels
            h=[];
            for i=1:height(link)
                if iscell(link.source(i))
                    source = link.source{i};
                    detector = link.detector{i};
                else
                    source = link.source(i);
                    detector = link.detector(i);
                end
                for j=1:length(source)
                    s = source(j);
                    d = detector(j);

                    src = SrcPos3D(s,:);
                    det = DetPos3D(d,:);

                    [x,y,z] = cylinder2P(2,100,src,det);
                    surf(x,y,z,[],'FaceColor',[0 1 0],'EdgeAlpha',0);

                end
            end

            axis(axis_handle,'equal');

            if(nargout>0)
                varargout{1}=h;
            end

        end

        function varargout=draw1020interp(obj,colors, lineStyles, axis_handle)
            if isnumeric(obj.link.type)
                link = obj.link(obj.link.type==obj.link.type(1),1:2);
            else
                link = obj.link(strcmp(obj.link.type,obj.link.type(1)),1:2);
            end
            n = height(link);
            if nargin < 2 || isempty(colors)
                colors = repmat([0.3 0.5 1], [n 1]);
            elseif size(colors,1) == 1
                colors = repmat(colors, [n 1]);
            end

            if nargin < 3 || isempty(lineStyles)
                lineStyles = repmat({'LineStyle', '-', 'LineWidth', 6}, [n 1]);
            elseif size(lineStyles, 1) == 1
                lineStyles = repmat({'LineStyle', '-', 'LineWidth', 6}, [n 1]);
            end

            if nargin < 4
                axis_handle = axes();
            end

            h=draw1020(obj,[],[],axis_handle);
            delete(h);
            lstS=find(ismember(obj.optodes.Type,'Source'));
            lstD=find(ismember(obj.optodes.Type,'Detector'));

            lst=find(~ismember(obj.optodes_registered.Type,{'FID-anchor','FID-attractor'}));
            Pos(:,1)=obj.optodes_registered.X(lst);
            Pos(:,2)=obj.optodes_registered.Y(lst);
            Pos(:,3)=obj.optodes_registered.Z(lst);

            [x,y]=obj.convert2d(Pos);
            xylink = [];
            xycolors = [];
            for i=1:height(link)
                if iscell(link.source(i))
                    source = link.source{i};
                    detector = link.detector{i};
                else
                    source = link.source(i);
                    detector = link.detector(i);
                end
                for j = 1:length(source)
                    s=source(j);
                    d=detector(j);
                    xylink(end+1,1)=mean(x([lstS(s) lstD(d)]));
                    xylink(end,2)=mean(y([lstS(s) lstD(d)]));
                    xycolors(end+1,:) = colors(i,:);
                end
            end

            xlim=get(axis_handle,'XLim');
            ylim=get(axis_handle,'YLim');
            [x2,y2]=meshgrid(xlim(1):xlim(2),ylim(1):ylim(2));

            F = scatteredInterpolant(xylink(:,1),xylink(:,2),xycolors(:,1),'nearest','none');
            cm(:,:,1)=reshape(F(x2(:),y2(:)),size(x2));
            F.Values=xycolors(:,2);
            cm(:,:,2)=reshape(F(x2(:),y2(:)),size(x2));
            F.Values=xycolors(:,3);
            cm(:,:,3)=reshape(F(x2(:),y2(:)),size(x2));


            k=dsearchn(xycolors,reshape(cm,[],3));
            k=reshape(k,size(cm,1),size(cm,2));
            lst=find(ismember(k,find(strcmp(lineStyles(:,2),'--'))));

            mask=ones(size(cm,1),size(cm,2));
            mask(lst)=NaN;
            cm(:,:,1)=mask.*cm(:,:,1);
            cm(:,:,2)=mask.*cm(:,:,2);
            cm(:,:,3)=mask.*cm(:,:,3);

            i=imagesc(axis_handle,x2(:),y2(:),cm);
            h=draw1020(obj,[],[],axis_handle);
            set(h,'LineWidth',.1,'color',[.3 .3 .3]);
            set(i,'alphaData',~isnan(cm(:,:,1)));

            if(obj.zoom)
                dx=(max(x)-min(x))/20;
                dy=(max(y)-min(y))/20;
                set(gca,'Xlim',[min(x)-dx max(x)+dx]);
                set(gca,'Ylim',[min(y)-dy max(y)+dy]);
            end

            if(nargout>0)
                varargout{1}=h;
            end
        end

        function varargout=draw1020(obj,colors, lineStyles, axis_handle,addlabels)
            % Code to draw the probe in 10-20 space

            if(~isempty(obj.link))
                if isnumeric(obj.link.type)
                    link = obj.link(obj.link.type==obj.link.type(1),1:2);
                else
                    link = obj.link(strcmp(obj.link.type,obj.link.type(1)),1:2);
                end
            else
                link=table;
            end

            if(isa(link,'table'))
                n = height(link);
            else
                n=size(link,1);
            end

            if nargin < 2 || isempty(colors)
                colors = repmat([0.3 0.5 1], [n 1]);
            elseif size(colors,1) == 1
                colors = repmat(colors, [n 1]);
            end

            if nargin < 3 || isempty(lineStyles)
                lineStyles = repmat({'LineStyle', '-', 'LineWidth', 3}, [n 1]);
            elseif size(lineStyles, 1) == 1
                lineStyles = repmat({'LineStyle', '-', 'LineWidth', 3}, [n 1]);
            end

            if nargin < 4 || isempty(axis_handle)
                axis_handle = axes();
            end
            if(nargin<5)
                addlabels=false;
            end

            hold(axis_handle,'on');
            set(axis_handle,'Projection','orthographic');
            view(axis_handle,[0 -90]);

            [x,y]=obj.convert2d(obj.pts1020);
            dx=-x(find(ismember(obj.labels,'Cz')));
            dy=-y(find(ismember(obj.labels,'Cz')));
            scatter(x+dx,y+dy,'filled','MarkerFaceColor',[.8 .8 .8],'parent',axis_handle);

            %             tbl=nirs.util.list_1020pts('?');
            %             for i=1:height(tbl); t(i)=text(x(i)+dx,y(i)+dy,tbl.Name{i}); end;
            %

            % Todo-  draw the probe too
            h=[];
            if(~isempty(obj.optodes_registered))
                % Points from the probe
                lst=find(~ismember(obj.optodes_registered.Type,{'FID-anchor','FID-attractor'}));
                Pos(:,1)=obj.optodes_registered.X(lst);
                Pos(:,2)=obj.optodes_registered.Y(lst);
                Pos(:,3)=obj.optodes_registered.Z(lst);
                [x,y]=obj.convert2d(Pos);

                DetPos3D=obj.detPos3D;
                SrcPos3D=obj.srcPos3D;
                [Detx,Dety]=obj.convert2d(DetPos3D);
                [Srcx,Srcy]=obj.convert2d(SrcPos3D);

                xop=x; yop=y;
                lstS=find(ismember(obj.optodes_registered.Type,'Source'));
                scatter(x(lstS)+dx,y(lstS)+dy,'filled','MarkerFaceColor','r','parent',axis_handle)
                lstD=find(ismember(obj.optodes_registered.Type,'Detector'));
                scatter(x(lstD)+dx,y(lstD)+dy,'filled','MarkerFaceColor','b','parent',axis_handle)

                if(addlabels)
                    for i=1:length(lstS)
                        xx=x(lstS(i))+dx;
                        yy=y(lstS(i))+dy;
                        text(xx, yy,['S' num2str(i)],'HorizontalAlignment','center','parent',axis_handle, 'FontSize', 12);
                    end
                    for i=1:length(lstD)
                        xx=x(lstD(i))+dx;
                        yy=y(lstD(i))+dy;
                        text(xx, yy,['D' num2str(i)],'HorizontalAlignment','center','parent',axis_handle, 'FontSize', 12);
                    end
                end

                for i=1:height(link)
                    if iscell(link.source(i))
                        source = link.source{i};
                        detector = link.detector{i};
                    else
                        source = link.source(i);
                        detector = link.detector(i);
                    end
                    for j=1:length(source)
                        s = source(j);
                        d = detector(j);
                        h(i)=line([Srcx(s) Detx(d)]+dx,[Srcy(s) Dety(d)]+dy,'parent',axis_handle,'Color', colors(i, :), lineStyles{i, :});
                        set(h(i),'UserData',[s d]);
                    end
                end
                %%
            else
                h=[];
            end


            [x,y]=obj.convert2d(obj.pts1020);
            headradius=obj.headcircum/(2*pi);

            headradius=norm([x(find(ismember(obj.labels,'Fpz'))) y(find(ismember(obj.labels,'Fpz')))]-...
                [x(find(ismember(obj.labels,'Cz'))) y(find(ismember(obj.labels,'Cz')))]);

            %              headradius=norm(obj.pts1020(find(ismember(obj.labels,'Fpz')),1:2)-...
            %                  obj.pts1020(find(ismember(obj.labels,'Cz')),1:2));
            %
            % add a circle for the head
            theta = linspace(0,2*pi);
            plot(headradius*cos(theta),headradius*sin(theta),'color',[.4 .4 .4],'linestyle','--','parent',axis_handle);

            headradius=norm([x(find(ismember(obj.labels,'nas'))) y(find(ismember(obj.labels,'nas')))]-...
                [x(find(ismember(obj.labels,'Cz'))) y(find(ismember(obj.labels,'Cz')))]);

            plot(headradius*cos(theta),headradius*sin(theta),'k','parent',axis_handle);
            plot([-headradius headradius],[0 0],'color',[.6 .6 .6],'linestyle','--','parent',axis_handle)
            plot([0 0],[-headradius headradius],'color',[.6 .6 .6],'linestyle','--','parent',axis_handle)

            line([-10 0],[-headradius -headradius-10],'color','k','parent',axis_handle);
            line([10 0],[-headradius -headradius-10],'color','k','parent',axis_handle);
            scatter([-15 15],[-headradius -headradius],'filled','k','sizedata',120,'parent',axis_handle);

            % Draw the central sulcus

            %lst={'CpZ','C3','C5-FC5'}
            pts=[x(find(ismember(obj.labels,'CPz'))), x(find(ismember(obj.labels,'C3'))), ...
                (x(find(ismember(obj.labels,'C5')))+x(find(ismember(obj.labels,'FC5'))))/2;...
                y(find(ismember(obj.labels,'CPz'))), y(find(ismember(obj.labels,'C3'))), ...
                (y(find(ismember(obj.labels,'C5')))+y(find(ismember(obj.labels,'FC5'))))/2]';
            pts(:,1)=pts(:,1)+dx; pts(:,2)=pts(:,2)+dy;
            xx=[min(pts(:,1)):.1:max(pts(:,1))];
            p=polyfit(pts(:,1),pts(:,2),2);
            plot(xx,polyval(p,xx),'color',[.4 .4 .4],'linestyle','-','parent',axis_handle);

            %lst2={'CpZ','C4','C6-FC6'}
            pts=[x(find(ismember(obj.labels,'CPz'))), x(find(ismember(obj.labels,'C4'))), ...
                (x(find(ismember(obj.labels,'C6')))+x(find(ismember(obj.labels,'FC6'))))/2;...
                y(find(ismember(obj.labels,'CPz'))), y(find(ismember(obj.labels,'C4'))), ...
                (y(find(ismember(obj.labels,'C6')))+y(find(ismember(obj.labels,'FC6'))))/2]';
            pts(:,1)=pts(:,1)+dx; pts(:,2)=pts(:,2)+dy;
            xx=[min(pts(:,1)):.1:max(pts(:,1))];
            p=polyfit(pts(:,1),pts(:,2),2);
            plot(xx,polyval(p,xx),'color',[.4 .4 .4],'linestyle','-','parent',axis_handle);

            % Add the insular sulcus
            % lst={'FT9','FT7','C5','Cp5'}

            pts=[x(find(ismember(obj.labels,'FT9'))), x(find(ismember(obj.labels,'FT7'))), ...
                x(find(ismember(obj.labels,'C5'))), x(find(ismember(obj.labels,'CP5')));...
                y(find(ismember(obj.labels,'FT9'))), y(find(ismember(obj.labels,'FT7'))), ...
                y(find(ismember(obj.labels,'C5'))), y(find(ismember(obj.labels,'CP5')))]';
            pts(:,1)=pts(:,1)+dx; pts(:,2)=pts(:,2)+dy;
            yy=[min(pts(:,2)):.1:max(pts(:,2))];
            p=polyfit(pts(:,2),pts(:,1),3);
            plot(polyval(p,yy),yy,'color',[.4 .4 .4],'linestyle','-','parent',axis_handle);

            % lst={'FT10','FT8','C6','Cp6'}
            pts=[x(find(ismember(obj.labels,'FT10'))), x(find(ismember(obj.labels,'FT8'))), ...
                x(find(ismember(obj.labels,'C6'))), x(find(ismember(obj.labels,'CP6')));...
                y(find(ismember(obj.labels,'FT10'))), y(find(ismember(obj.labels,'FT8'))), ...
                y(find(ismember(obj.labels,'C6'))), y(find(ismember(obj.labels,'CP6')))]';
            pts(:,1)=pts(:,1)+dx; pts(:,2)=pts(:,2)+dy;
            yy=[min(pts(:,2)):.1:max(pts(:,2))];
            p=polyfit(pts(:,2),pts(:,1),3);
            plot(polyval(p,yy),yy,'color',[.4 .4 .4],'linestyle','-','parent',axis_handle);

            %set(gcf,'color','w');

            if(obj.zoom)
                dx=(max(xop)-min(xop))/10;
                dy=(max(yop)-min(yop))/10;
                set(axis_handle,'Xlim',[min(xop)-dx max(xop)+dx]);
                set(axis_handle,'Ylim',[-headradius*1.13 max(yop)+dy]);
            end
            set(axis_handle,'Xdir','reverse');
            %set(axis_handle,'YDir','reverse')
            if(nargout>0)
                varargout{1}=h;
            end
            axis(axis_handle,'tight');
            axis(axis_handle,'equal');
            axis(axis_handle,'off');
            set(get(axis_handle,'parent'),'color','w');
        end

        function headcircum = get.headcircum(obj)
            a=obj.LR_distance/2;
            b=obj.AP_distance/2;
            headcircum=.9*pi*(3*(a+b)-sqrt((3*a+b)*(a+3*b)));

        end

        function LR_arclength = get.LR_arclength(obj)
            a=obj.LR_distance/2;
            b=obj.IS_distance;
            LR_arclength=pi*(3*(a+b)-sqrt((3*a+b)*(a+3*b)))/2;

        end

        function AP_arclength = get.AP_arclength(obj)
            a=obj.AP_distance/2;
            b=obj.IS_distance;
            AP_arclength=pi*(3*(a+b)-sqrt((3*a+b)*(a+3*b)))/2;

        end

        function [x,y]=convert2d(obj,pts)

            %Clarke far-side general prospective azumuthal projection
            r = -2.4;
            R=sqrt(sum(pts.^2,2));
            x=r*R.*(pts(:,1)./abs(pts(:,3)-r*R));
            y=r*R.*(pts(:,2)./abs(pts(:,3)-r*R));

        end

        function fwdModel = convert2FwdModel(obj,type)
            if(nargin<2)
                type='SLAB';
            end
            switch(upper(type))
                case('BEM')
                    fwdModel=nirs.forward.NirfastBEM;
                    fwdModel.mesh=obj.getmesh;
                    fwdModel.prop=obj.opticalproperties;
                    fwdModel.probe=obj;
                case('FEM')
                    fwdModel=nirs.forward.NirfastFEM;
                    fwdModel.mesh=obj.getmesh;
                    fwdModel.mesh=fwdModel.mesh.convert2FEM;
                    fwdModel.prop=obj.opticalproperties;
                    fwdModel.probe=obj;
                case('SLAB')
                    fwdModel=nirs.forward.ApproxSlab;
                    fwdModel.mesh=obj.getmesh;
                    fwdModel.prop=obj.opticalproperties;
                    fwdModel.probe=obj;
                case('MMC')
                    fwdModel=nirs.nirs.forward.MMCLab;
                    fwdModel.mesh=obj.getmesh;
                    fwdModel.mesh=fwdModel.mesh.convert2FEM;
                    fwdModel.prop=obj.opticalproperties;
                    fwdModel.probe=obj;
                case('MCX')
                    fwdModel=nirs.forward.MCXLab;
                    fwdModel.mesh=obj.getmesh;
                    fwdModel.mesh=fwdModel.mesh.convert2image;
                    fwdModel.prop=obj.opticalproperties;
                    fwdModel.probe=obj;
            end

        end


        %draw1020image;  % Code to draw an image in 10-20space
        % makeimage;  % code to do simple image reconstruction in 10-20 space
        %registerprobe;
        % depth map
        % define ROIs
        % coregister
        % register probe



    end

end