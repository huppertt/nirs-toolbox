classdef ProbeHyperscan1020 < nirs.core.ProbeHyperscan
    %% ProbeHyperscan.  This is a class derived from the core Probe class which adds two-person drawing

    %% PROBE - This object hold nirs probe geometries.
    %
    % Properties:
    %     optodes-  a table containing the position, name, and type of any
    %     probe points including any reference or anchor points used in
    %     registration
    %     srcPos - (dependent) nsrc x 3 array of source positions (mm)
    %     detPos - (dependent) ndet x 3 array of detector positions (mm)
    %     link   - a table containing the columns 'source', 'detector',
    %              and 'type' describing the the connections between sources
    %              and detectors
    %     distances - (dependent) returns measurement distances
    %
    %  Methods:
    %     swapSD - returns a Probe object with sources exchanged for detectors
    %     draw   - displays a visualization of the probe geometry

    properties
        defaultdrawfcn='10-20'; % Drawing function to use
    end

    properties( Dependent = true )
        optodes_registered;   % Registered version of optodes
        srcPos3D;      % nsrc x 3 array of source positions (mm)
        detPos3D;      % ndet x 3 array of detector positions (mm)
        headcircum;  % circumference of head  % make dependent
        AP_arclength; % Anterior-Posterior arc length
        LR_arclength; % Left-right arclength
        braindepth;  % depth of brain for modeling
        opticalproperties;  % optical properties for modeling
          srcPos_drawing;      % nsrc x 3 array of source positions (mm)
        detPos_drawing;      % ndet x 3 array of detector positions (mm)
    end

    properties(Hidden = true)
        RotateMatrix3D;
        zoom; % Flag to zoom in or show full 10-20 probe

        AP_distance;  % Distance from Oz to Fpz
        LR_distance;  % Distance from LPA to RPA
        IS_distance;  % Distance from center to Cz
        labels;  % labels of 10-20 points
        pts1020;  % 10-20 points
        mesh;
      
    end

    methods
      
        function obj = set_mesh(obj,mesh)
            if(~iscell(mesh))
                mesh={mesh};
            end
            if(length(mesh)~=length(obj.originalprobe))
                mesh=repmat(mesh,length(obj.originalprobe),1);
            end

            for i=1:length(obj.originalprobe)
                obj.originalprobe(i)=obj.originalprobe(i).set_mesh(mesh{i});
            end
        end
        function obj = SetFiducialsVisibility(obj,flag)
            for i=1:length(obj.originalprobe)
                obj.originalprobe(i)=obj.originalprobe(i).SetFiducialsVisibility(flag);
            end
        end

        function mesh = getmesh(obj)
             for i=1:length(obj.originalprobe)
                 mesh{i}=obj.originalprobe(i).getmesh;
             end
        end


        
        function detPos_drawing = get.detPos_drawing(obj)
            detPos_drawing=[];
            c=get(gca,'children');
            for i=1:length(obj.SubjectLabels)
                for j=1:length(c)
                    if(strcmp(get(c(j),'Tag'),[obj.SubjectLabels{i} '-DetectorOptodes']))
                        if(isfield(c(j),'ZData'))
                            detPos_drawing= [detPos_drawing; [get(c(j),'XData'); get(c(j),'YData'); get(c(j),'ZData')]'];

                        else
                            detPos_drawing= [detPos_drawing; [get(c(j),'XData'); get(c(j),'YData'); zeros(size(get(c(j),'YData')))]'];
                        end
                    end
                end
            
            end
            if(isempty(detPos_drawing))
                % If it was drawn from a nirs.core.Probe, there will
                % not be SubjectLabels info
                for j=1:length(c)
                    if(strcmp(get(c(j),'Tag'),'DetectorOptodes'))
                        if(isfield(c(j),'ZData'))
                            detPos_drawing= [detPos_drawing; [get(c(j),'XData'); get(c(j),'YData'); get(c(j),'ZData')]'];

                        else
                            detPos_drawing= [detPos_drawing; [get(c(j),'XData'); get(c(j),'YData'); zeros(size(get(c(j),'YData')))]'];
                        end
                    end
                end
            end
        end
        function srcPos_drawing = get.srcPos_drawing(obj)
            srcPos_drawing=[];
            c=get(gca,'children');
            for i=1:length(obj.SubjectLabels)
                for j=1:length(c)
                    if(strcmp(get(c(j),'Tag'),[obj.SubjectLabels{i} '-SourceOptodes']))
                           if(isfield(c(j),'ZData'))
                            srcPos_drawing= [srcPos_drawing; [get(c(j),'XData'); get(c(j),'YData'); get(c(j),'ZData')]'];

                        else
                            srcPos_drawing= [srcPos_drawing; [get(c(j),'XData'); get(c(j),'YData'); zeros(size(get(c(j),'YData')))]'];
                        end
                    end
                end
                
            end
            if(isempty(srcPos_drawing))
                % If it was drawn from a nirs.core.Probe, there will
                % not be SubjectLabels info
                for j=1:length(c)
                    if(strcmp(get(c(j),'Tag'),'SourceOptodes'))
                        if(isfield(c(j),'ZData'))
                            srcPos_drawing= [srcPos_drawing; [get(c(j),'XData'); get(c(j),'YData'); get(c(j),'ZData')]'];

                        else
                            srcPos_drawing= [srcPos_drawing; [get(c(j),'XData'); get(c(j),'YData'); zeros(size(get(c(j),'YData')))]'];
                        end
                    end
                end
            end
        end
        function obj=set.defaultdrawfcn(obj,defaultdrawfcn)
            obj.defaultdrawfcn=defaultdrawfcn;
            for i=1:length(obj.originalprobe)
                obj.originalprobe(i).defaultdrawfcn=obj.defaultdrawfcn;
            end
        end


        function braindepth = get.braindepth(obj)
            for i=1:length(obj.originalprobe)
                braindepth(i)=obj.originalprobe(i).braindepth;
            end
        end
        function opticalproperties = get.opticalproperties(obj)
            for i=1:length(obj.originalprobe)
                opticalproperties(i)=obj.originalprobe(i).opticalproperties;
            end
        end

        function srcPos3D = get.srcPos3D(obj)
            %% This function returns the src pos (in mm)
            [tbl,l]=sortrows(obj.optodes_registered,{'Type','Name'});
            lst=find(ismember(tbl.Type,'Source'));
            srcPos3D=[tbl.X(lst) tbl.Y(lst) tbl.Z(lst)];

            %Convert to mm if needed
            lstCM=find(ismember(tbl.Units(lst),{'cm'}));
            srcPos3D(lstCM,:)=srcPos3D(lstCM,:)*10;
        end

       
        function optodes_registered = get.optodes_registered(obj)
            optodes_registered=table;
            curSrc=0;
            curDet=0;
            for i=1:length(obj.originalprobe)
                opttmp=obj.originalprobe(i).optodes_registered;
                xyz=[opttmp.X opttmp.Y opttmp.Z];
                xyz(:,4)=1;
                xyz=xyz*obj.RotateMatrix{i};
                opttmp.X=xyz(:,1);
                opttmp.Y=xyz(:,2);
                opttmp.Z=xyz(:,3);
                opttmp.SubjectLabel=repmat(obj.SubjectLabels{i},height(opttmp),1);

                ltmp=obj.originalprobe(i).link;
                ltmp.source=ltmp.source+curSrc;
                ltmp.detector=ltmp.detector+curDet;


                for j=1:height(opttmp)
                    if(ismember(opttmp.Type{j},{'Source'}))
                        n=opttmp.Name{j};
                        n=str2num(strrep(n,'Source-',''));
                        n=['0000' num2str(n+curSrc)];
                        n=['Source-' n(end-3:end)];
                        opttmp.Name{j}=n;

                    elseif(ismember(opttmp.Type{j},{'Detector'}))
                        n=opttmp.Name{j};
                        n=str2num(strrep(n,'Detector-',''));
                        n=['0000' num2str(n+curDet)];
                        n=['Detector-' n(end-3:end)];
                        opttmp.Name{j}=n;
                    end

                end
                curDet=max(ltmp.detector);
                curSrc=max(ltmp.source);

                optodes_registered=[optodes_registered; opttmp];
            end

        end

        function varargout=draw( obj, varargin)
             for i=1:length(obj.originalprobe)
                obj.originalprobe(i).defaultdrawfcn=obj.defaultdrawfcn;
             end
            
             

            if(strcmp(lower(obj.defaultdrawfcn),'2d'))
                optodes= obj.optodes;
                h={};
                for i=1:length(obj.SubjectLabels)
                    pp=nirs.core.Probe;
                    pp.link=obj.originalprobe(i).link;
                    pp.optodes=optodes(ismember(optodes.SubjectLabel,obj.SubjectLabels{i}),:);
                    h{i}=pp.draw([],[],axis_handle,true);
                end
                
                c=get(gca,'children');
              


            else
                hold on;
                set(get(gca,'children'),'Tag','none');
                if(length(varargin)<3 || ~isa(varargin{3},'matlab.graphics.axis.Axes'))
                    varargin{3}=gca;
                end
                for i=1:length(obj.originalprobe)
                    hl{i}=obj.originalprobe(i).draw(varargin{:});
                    
                    if(obj.show_labels)
                        % make sure text is always above the object
                        ang=atan2(norm(cross([1 0 0],[0 1 0]*obj.RotateMatrix{i}(1:3,1:3))),dot([1 0 0],[0 1 0]*obj.RotateMatrix{1}(1:3,1:3)));

                        if(contains(lower(obj.defaultdrawfcn),'3d'))
                            com=[mean(get(varargin{3},'Xlim')) mean(get(varargin{3},'Ylim')) 1.5*max(get(varargin{3},'Zlim'))];
                        else
                            com=[mean(get(varargin{3},'Xlim')) 1.5*max(get(varargin{3},'Ylim')) 0];
                        end
                        text(varargin{3},com(1),com(2),com(3),obj.SubjectLabels{i},...
                            'FontWeight','bold','FontSize',16,'Rotation',ang/pi*180,...
                            'HorizontalAlignment','center');
                    end

                    % rescale(i,1)=mean(obj.originalprobe(i).optodes.X);
                    % rescale(i,2)=mean(obj.originalprobe(i).optodes.Y);
                    c=get(varargin{3},'children');

                    for j=1:length(c)
                        if(isempty(get(c(j),'Tag')))
                            set(c(j),'Tag',obj.SubjectLabels{i});
                        elseif(strcmp(get(c(j),'Tag'),'DetectorOptodes'))
                            set(c(j),'Tag',[obj.SubjectLabels{i} '-DetectorOptodes']);
                         elseif(strcmp(get(c(j),'Tag'),'SourceOptodes'))
                            set(c(j),'Tag',[obj.SubjectLabels{i} '-SourceOptodes']);
                        end
                    end
                end

                c=get(varargin{3},'children');

                for i=1:length(obj.originalprobe)
                    xyz=[];
                    for j=1:length(c)
                        if(strcmp(get(c(j),'Tag'),obj.SubjectLabels{i}) |...
                                strcmp(get(c(j),'Tag'),[obj.SubjectLabels{i} '-DetectorOptodes']) | ...
                                strcmp(get(c(j),'Tag'),[obj.SubjectLabels{i} '-SourceOptodes']))

                            if(isa(c(j),'matlab.graphics.primitive.Text'))
                                xyz=[xyz; get(c(j),'Position')];
                            elseif(isa(c(j),'matlab.graphics.primitive.Patch'))
                                xyz=[xyz; get(c(j),'Vertices')];


                            elseif(isa(c(j),'matlab.graphics.chart.primitive.Scatter') | ...
                                    isa(c(j),'matlab.graphics.primitive.Line') | ...
                                    isa(c(j),'matlab.graphics.chart.primitive.Line'))

                                x=get(c(j),'XData');
                                y=get(c(j),'YData');

                                z=get(c(j),'ZData');
                                if(isempty(z))
                                    xyz2=[x(:) y(:)];
                                    xyz2(:,3)=0;
                                else
                                    xyz2=[x(:) y(:) z(:)];
                                end
                                xyz=[xyz; xyz2];

                            end
                        end
                    end
                    recenter(i,:)=(max(xyz,[],1)+min(xyz,[],1))/2;
                    scale(i)=sqrt(max(sum((xyz-ones(size(xyz,1),1)*recenter(i,:)).^2,2)))*1.5;
                end


                c=get(varargin{3},'children');
                for i=1:length(obj.originalprobe)
                    Rot=obj.RotateMatrix{i};
                    if(contains(lower(obj.defaultdrawfcn),'3d'))

                        % The 3D axis directions are flipped
                        a=pi;
                        T=[cos(a) -sin(a) 0 0; sin(a) cos(a) 0 0; 0 0 1 0; 0 0 0 1];
                        Rot=T*Rot;

                        Rot(4,1)=max(min(Rot(4,1),2),-2)*scale(i);
                        Rot(4,2)=max(min(Rot(4,2),2),-2)*scale(i);
                        Rot(4,3)=max(min(Rot(4,3),2),-2)*scale(i);
                    else
                        Rot(4,1)=max(min(Rot(4,1),1),-1)*scale(i);
                        Rot(4,2)=max(min(Rot(4,2),1),-1)*scale(i);
                        Rot(4,3)=max(min(Rot(4,3),1),-1)*scale(i);
                    end
                    for j=1:length(c)

                        if(strcmp(get(c(j),'Tag'),obj.SubjectLabels{i}) |...
                                strcmp(get(c(j),'Tag'),[obj.SubjectLabels{i} '-DetectorOptodes']) | ...
                                strcmp(get(c(j),'Tag'),[obj.SubjectLabels{i} '-SourceOptodes']))
                            
                            if(isa(c(j),'matlab.graphics.primitive.Text'))
                                xyz=get(c(j),'Position');
                                xyz(:,1)=xyz(:,1)-recenter(i,1);
                                xyz(:,2)=xyz(:,2)-recenter(i,2);
                                xyz(:,3)=xyz(:,3)-recenter(i,3);

                                xyz(:,4)=1;
                                xyz=xyz*Rot;
                                set(c(j),'Position',xyz(1:3));
                            elseif(isa(c(j),'matlab.graphics.primitive.Patch'))
                                v=get(c(j),'Vertices');
                                v(:,1)=v(:,1)-recenter(i,1);
                                v(:,2)=v(:,2)-recenter(i,2);
                                v(:,3)=v(:,3)-recenter(i,3);
                                v(:,4)=1;
                                v=v*Rot;
                                set(c(j),'Vertices',v(:,1:3));

                            elseif(isa(c(j),'matlab.graphics.chart.primitive.Scatter') | ...
                                    isa(c(j),'matlab.graphics.primitive.Line') | ...
                                    isa(c(j),'matlab.graphics.chart.primitive.Line'))

                                x=get(c(j),'XData'); 
                                y=get(c(j),'YData'); 

                                z=get(c(j),'ZData');
                                if(isempty(z))
                                    xyz=[x(:) y(:)];
                                    xyz(:,3)=0;
                                else
                                    xyz=[x(:) y(:) z(:)];
                                end
                                xyz(:,4)=1;
                                xyz(:,1)=xyz(:,1)-recenter(i,1);
                                xyz(:,2)=xyz(:,2)-recenter(i,2);
                                xyz(:,3)=xyz(:,3)-recenter(i,3);
                                
                                xyz=xyz*Rot;

                                set(c(j),'XData',xyz(:,1));
                                set(c(j),'YData',xyz(:,2));
                                if(~isempty(z))
                                    set(c(j),'ZData',xyz(:,3));
                                end
                            
                            end
                            
                        end
                       
                    end
                    
                end


            end
            
            for i=1:length(obj.SubjectLabels)
                pp.link=obj.originalprobe(i).link;
                ll{i}=pp.link(ismember(pp.link.type,pp.link.type(1)),:);
                x=zeros(length(hl{i}),1);
                y=zeros(length(hl{i}),1);
                z=zeros(length(hl{i}),1);

                for j=1:length(hl{i})
                    x(j)=mean(get(hl{i}(j),'Xdata'));
                    y(j)=mean(get(hl{i}(j),'Ydata'));
                    if(~isempty(get(hl{i}(j),'Zdata')))
                        z(j)=mean(get(hl{i}(j),'Zdata'));
                    else
                        z(j)=NaN;
                    end
                end
                pos{i}=[x y z];
            end

            link=obj.link;
            link=link(ismember(link.TypeOrigin,link.TypeOrigin(1)) & ismember(link.TypeDest,link.TypeOrigin(1)),:);

            hold(varargin{3},'on');
            hh=[];
            for i=1:height(link)
                p1=find(ismember(obj.SubjectLabels,link.SubjectLabelOrigin{i}));
                p2=find(ismember(obj.SubjectLabels,link.SubjectLabelDest{i}));
                l1=find(ll{p1}.source==link.sourceOrigin(i) & ...
                    ll{p1}.detector==link.detectorOrigin(i)); 
                l2=find(ll{p2}.source==link.sourceDest(i) & ...
                    ll{p2}.detector==link.detectorDest(i));
                xyz1=pos{p1}(l1,:);
                xyz2=pos{p2}(l2,:);
                if(isnan(xyz1(3)))
                    hh(i,1)=line([xyz1(1) xyz2(1)],[xyz1(2) xyz2(2)]);
                else
                    hh(i,1)=line([xyz1(1) xyz2(1)],[xyz1(2) xyz2(2)],[xyz1(3) xyz2(3)]);
                end
            end
            set(hh,'linewidth',1,'color',[.7 .7 .7]);
            hold(varargin{3},'off');
            
            if(nargout>0)
                varargout{1}=hh;
            end
        end
    end

end