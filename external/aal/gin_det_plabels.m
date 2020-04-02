function [Label2, Perc, nbpCluster2, PercentROI2, nbpROI2] = gin_det_plabels(XYZmm, MNIY, MNID, ROI,fact)
%
% Compute Labels and % for local maxima
% 
% 	XYZmm	: Local Maxima 
% 	MNID	: atlas volume with labels values
% 	ROI	: structure : link between labels values and label name.
% 	BORDER_XYZ	:  structure contains the regions edges in mm
% 	BORDER_V	:  structure contains the labels values
%
%	Label2 		: Labels List found
%	Perc 		: Percent of cluster's voxels belonging to Label2
%	nbp_cluster2	: Number of voxels in the cluster
%	Percent_ROI2	: Percent of ROI's voxels (ROI of Label2) belonging to the cluster
%	nbp_roi2 	: Number of voxels in ROI of Label2
%_______________________________________________________________________
%
% gin_det_plabels.m				B Landeau 20/02/2002
% gin_det_plabels.m				B Landeau 19/04/2004
%_______________________________________________________________________

% first part
% transformation from XYZmm coord to MNID volume coord

	mniXYZ   = MNID.mat \ [XYZmm; ones(1, size(XYZmm, 2))];

	MNIV = 0;
	for i=1:size(mniXYZ,2),
	%if mniXYZ(1,i)>0 & mniXYZ(2,i)>0 & mniXYZ(3,i)>0
		MNIV(i) = MNIY(round(mniXYZ(1,i)),round(mniXYZ(2,i)),round(mniXYZ(3,i)));
	end

nbpTotalCluster = size(XYZmm,2); % ou size(MNIV,2)

% link labels values - labels
	nb_roi = length(ROI);
	for i=1:length(MNIV),
		Label(i).Nom =  'OUTSIDE';
		ID(i) = -1.0;
		nbpROI(i) = 0;
	end
	
	for i=1:nb_roi,
		
		tmp=find(MNIV==ROI(i).ID);
		if ~isempty(tmp)	
			tmp2 = find(MNIY==ROI(i).ID);
			for j=1:length(tmp),
				Label(tmp(j)).Nom = ROI(i).Nom_L;
				ID(tmp(j)) = ROI(i).ID;
				nbpROI(tmp(j)) = size(tmp2,1);

			end
		end
	end

% sort labels
% size(y,2) = nbpTotalCluster
% idx = nb de voxels de l intersection cluster-roi

	[y,is] 		= sort(ID);
	Label1 		= Label(is);
	nbpROI1 	= nbpROI(is);

	[u,idx]		= unique(y);
	Label2		= Label1(idx);
	nbpROI2 	= nbpROI1(idx);

% sort %
	Perc(1) = (idx(1))/size(y,2)*100.0;
	for jj=2:size(u,2)
		Perc(jj)=(idx(jj)-idx(jj-1))/size(y,2)*100.0;
	end

	PercentROI(1) =0;
	if nbpROI2(1) ~= 0
	   PercentROI(1) = (idx(1))/(nbpROI2(1)*fact)*100.0;
	end
	nbpTCluster(1) = nbpTotalCluster;

	for jj=2:size(u,2)
	   PercentROI(jj) =0;
	   if nbpROI2(jj) ~= 0
		PercentROI(jj)=(idx(jj)-idx(jj-1))/(nbpROI2(jj)*fact)*100.0;
	   end
	   nbpTCluster(jj) = nbpTotalCluster;
	end

	

	[Perc,inx] 	= sort(Perc);
	Label2 		= Label2(inx);
	PercentROI2 	= PercentROI(inx);
	nbpCluster2 	= nbpTCluster(inx);
	nbpROI2 	= nbpROI2(inx);

