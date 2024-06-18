classdef ProbeHyperscan
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
        SubjectLabels;
        show_labels=false;;
    end

    properties( Dependent = true )
        optodes    % table describing the src/det and any additional probe points
        link        % table describing the connections of source/detector pairs

        distances   % (dependent) returns measurement distances
        srcPos      % nsrc x 3 array of source positions (mm)
        detPos      % ndet x 3 array of detector positions (mm)
    end

    properties(Hidden = true)
        fixeddistances;
        originalprobe;
        RotateMatrix;

    end

    methods

        function obj = ProbeHyperscan(probe,SubjectLabels,RotateMatrix)


            if(nargin>1)
                obj.SubjectLabels=SubjectLabels;
            else
                for i=1:length(probe)
                    obj.SubjectLabels{i,1}=char(64+i);
                end
            end
            if(length(probe)<length(obj.SubjectLabels))
                probe=repmat(probe,length(obj.SubjectLabels),1);
            end
            obj.originalprobe=probe;

            if(nargin>2)
                obj.RotateMatrix=RotateMartix;
            else
                obj.RotateMatrix=cell(length(unique(obj.SubjectLabels)),1);
                maxY=-inf;
                minY=inf;
                for i=1:length(probe);
                    maxY=max(maxY,max(probe(i).optodes.Y));
                    minY=min(minY,min(probe(i).optodes.Y));
                    maxY=max(maxY,max(probe(i).optodes.X));
                    minY=min(minY,min(probe(i).optodes.X));
                end
                y=(maxY-minY);

                for i=1:length(obj.RotateMatrix)
                    a=(i-1)*2*pi/length(obj.RotateMatrix);
                    T=eye(4);
                    T(4,1)=sin(a);
                    T(4,2)=cos(a);
                    obj.RotateMatrix{i}=[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]'*[cos(a) -sin(a) 0 0; sin(a) cos(a) 0 0; 0 0 1 0; 0 0 0 1]*T;




                end

            end

        end

        function obj=resize_images(obj,resize)
            % Because we are often dealing with parent-child, it is useful
            % to be able to resize the drawing without effecting the data

            % example:  probeHS=probeHS.resize_images([1 .8]);  makes the
            % 2nd subject 80% the size of the first

            if(length(resize)<length(obj.RotateMatrix))
                resize=repmat(resize,length(obj.RotateMatrix),1);
            end

            for i=1:length(obj.RotateMatrix)
                obj.RotateMatrix{i}(1:3,1:3)=obj.RotateMatrix{i}(1:3,1:3)*resize(i);
            end


        end

        function obj=setRotation(obj,rot,horizontal)
            % options are {'inward','toward',1) or ('outward','away',-1)

            if(nargin<3)
                horizontal=false;
            end

            if(isnumeric(rot))
                for i=1:length(rot)
                    rot2{i}=cellstr(num2str(rot(i)));
                end
                rot=rot2;
            end
            if(ischar(rot))
                rot=cellstr(rot);
            end

            if(length(rot)<length(obj.originalprobe))
                rot=repmat(rot,length(obj.originalprobe),1);
            end



            for i=1:length(rot)
                maxY=-inf;
                minY=inf;
                maxX=-inf;
                minX=inf;
                for ii=1:length(obj.originalprobe);
                    maxY=max(maxY,max(obj.originalprobe(ii).optodes.Y));
                    minY=min(minY,min(obj.originalprobe(ii).optodes.Y));
                    maxX=max(maxX,max(obj.originalprobe(ii).optodes.X));
                    minX=min(minX,min(obj.originalprobe(ii).optodes.X));
                end
                y=(maxY-minY);

                a=(i-1)*2*pi/length(obj.RotateMatrix);
                
                if(horizontal)
                    a=a+pi/2;
                end
                T=eye(4);
                T(4,1)=sin(a);
                T(4,2)=cos(a);
                
                
                
                if(ismember(lower(rot{i}),{'inward','toward','1'}))
                    a=a;
                elseif(ismember(lower(rot{i}),{'outward','away','-1'}))
                    a=a-pi;
                end

                obj.RotateMatrix{i}=[cos(a) -sin(a) 0 0; sin(a) cos(a) 0 0; 0 0 1 0; 0 0 0 1]*T;
            end


        end

        function srcPos = get.srcPos(obj)
            %% This function returns the src pos (in mm)
            [tbl,l]=sortrows(obj.optodes,{'Type','Name'});
            lst=find(ismember(tbl.Type,'Source'));
            srcPos=[tbl.X(lst) tbl.Y(lst) tbl.Z(lst)];

            %Convert to mm if needed
            lstCM=find(ismember(tbl.Units(lst),{'cm'}));
            srcPos(lstCM,:)=srcPos(lstCM,:)*10;
        end

        function detPos = get.detPos(obj)
            %% This function returns the det pos (in mm)
            [tbl,l]=sortrows(obj.optodes,{'Type','Name'});
            lst=find(ismember(tbl.Type,'Detector'));

            detPos=[tbl.X(lst) tbl.Y(lst) tbl.Z(lst)];
            %Convert to mm if needed
            lstCM=find(ismember(tbl.Units(lst),{'cm'}));
            detPos(lstCM,:)=detPos(lstCM,:)*10;
        end




        function optodes = get.optodes(obj)
            optodes=table;
            curSrc=0;
            curDet=0;
            for i=1:length(obj.originalprobe)
                opttmp=obj.originalprobe(i).optodes;
                xyz=[opttmp.X opttmp.Y opttmp.Z];
                xyz(:,4)=1;
                xyz(:,1)=xyz(:,1)-mean(xyz(:,1));
                xyz(:,2)=xyz(:,2)-mean(xyz(:,2));
                xyz(:,3)=xyz(:,3)-mean(xyz(:,3));
                
                scale=sqrt(max(sum(xyz(:,1:3).^2,2)))*.5;
                Rot=obj.RotateMatrix{i};
                Rot(4,1)=max(min(Rot(4,1),1),-1)*scale;
                        Rot(4,2)=max(min(Rot(4,2),1),-1)*scale;
                        Rot(4,3)=max(min(Rot(4,3),1),-1)*scale;

                xyz=xyz*Rot;
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

                optodes=[optodes; opttmp];
            end

        end


        function link =get.link(obj)
            link=table;
            curSrc=0;
            curDet=0;
            for i=1:length(obj.originalprobe)
                ltmp=obj.originalprobe(i).link;
                ltmp.source=ltmp.source+curSrc;
                ltmp.detector=ltmp.detector+curDet;
                ltmp.SubjectLabel=repmat(cellstr(obj.SubjectLabels{i}),height(ltmp),1);
                curDet=max(ltmp.detector);
                curSrc=max(ltmp.source);
                link=[link; ltmp];
            end

        end

        function varargout=draw( obj, varargin)
            if(nargout>0)
                varargout{1}=[];
            end

            pp=nirs.core.Probe;
            pp.link=obj.link;
            pp.optodes=obj.optodes;

            if(nargout==0)
                pp.draw(varargin{:});
            elseif(nargout>0)
                varargout{1}=pp.draw(varargin{:});
            end
            
        end



        function d = get.distances( obj )
            %% distances - Calculates measurement distance for each channel.

            if(isa(obj,'nirs.core.ProbeROI'))
                d=obj.fixeddistances;
                return
            end


            isrc = obj.link.source;
            idet = obj.link.detector;

            if iscell(isrc)
                for i = 1:length(isrc)
                    vec = obj.srcPos(isrc{i},:) - obj.detPos(idet{i},:);
                    d(i,1) = mean(sqrt(sum(vec.^2,2)));
                end
            else
                vec = obj.srcPos(isrc,:) - obj.detPos(idet,:);
                d = sqrt( sum( vec.^2,2 ) );
            end

            if(~isempty(obj.fixeddistances))
                if(any((d-obj.fixeddistances)~=0))
                    warning('probe distances do not match layout');
                end
                d=obj.fixeddistances;
                return
            end


        end

        function obj = swapSD( obj )
            %% swapSD - Swaps sources for detectors and vice versa.

            lstD=find(ismember(obj.optodes.Type,'Detector'));
            lstS=find(ismember(obj.optodes.Type,'Source'));
            obj.optodes.Type(lstS)=repmat({'Detector'},length(lstS),1);
            obj.optodes.Type(lstD)=repmat({'Source'},length(lstD),1);


            obj.link(:,[1 2]) = obj.link(:, [2 1]);
        end


    end

end