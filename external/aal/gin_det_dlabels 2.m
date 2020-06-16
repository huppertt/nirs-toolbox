function [Label,Dist] = gin_det_dlabels(XYZmm, MNID, ROI, BORDER_XYZ, BORDER_V)
%
% Compute Labels and Distances for local maxima
% 
% 	XYZmm	: Local Maxima 
% 	MNID	: atlas volume with labels values
% 	ROI	: structure : link between labels values and label name.
% 	BORDER_XYZ	:  structure contains the regions edges in mm
% 	BORDER_V	:  structure contains the labels values
%
%	Label : Labels List found
%	Dist  : Distances List for Labels List
%
%_______________________________________________________________________
%
% gin_det_plabels.m				B Landeau 20/02/02
%_______________________________________________________________________


% first part
% Determine the Labels & distances(0.0) for local maxima in region
% transformation from XYZmm coord to MNID volume coord
 	mniXYZ   = MNID.mat \ [XYZmm; ones(1, size(XYZmm, 2))];
 	MNIY     = spm_read_vols(MNID);
 	for i=1:size(mniXYZ,2),
 		%MNIV(i) = MNIY(mniXYZ(1,i),mniXYZ(2,i),mniXYZ(3,i));
		MNIV(i) = MNIY(round(mniXYZ(1,i)),round(mniXYZ(2,i)),round(mniXYZ(3,i)));
 	end

% link labels values - labels
 	nb_roi = length(ROI);
 	for i=1:length(MNIV),
 		Label(i).Nom = 'OUTSIDE';
 		Dist(i) = -1.0;
		dd=Dist(i);	% <- ajout le 18/07/01
 	end
 	for i=1:nb_roi,
 		tmp=find(MNIV==ROI(i).ID);
 		if ~isempty(tmp)	
 			for j=1:length(tmp),
 				Label(tmp(j)).Nom = ROI(i).Nom_L;
 				Dist(tmp(j)) = 0.0;
				nn=Label(tmp(j)).Nom;	% <- ajout le 18/07/01
				dd=Dist(tmp(j));		% <- ajout le 18/07/01
 			end
 		end
 	end
	

% part two
% Determine the Labels & Distances for local maxima outside regions
% if Dist(1)==-1.0   <- modification le 18/07/01

	% transformation from XYZmm coord to MNID volume coord
	xyzmm = [XYZmm;1] * ones(1, size(BORDER_XYZ, 2));
	bxyzmm = MNID.mat * [BORDER_XYZ; ones(1, size(BORDER_XYZ, 2))];

	% calculate distance 
	d=sqrt(sum((xyzmm-bxyzmm).^2,1));

	% sort distance 
	[ds,index]=sort(d);
	val=BORDER_V(index);

	value(1)=0;
	value(2)=0;
	value(3)=0;
	ivalue(1)=0;
	ivalue(2)=0;
	ivalue(3)=0;

	jj=1;
	for k=1:size(val,2),
		if jj==4
			break;
		else
			for j=1:jj,
				if val(k)==value(j)
					found=1;
					break;
				else
					found=0;
				end
			end
			if found==0
				value(jj)=val(k);
				ivalue(jj)=k;
				jj=jj+1;
			end
		end
	end

% keep the smallest 3 distances
	Dist=ds(ivalue);

% search the link values - labels
	nb_roi = length(ROI);

	for i=1:length(value),
		Label(i).Nom = 'OUTSIDE';
	end

	for i=1:nb_roi,
		tmp=find(value==ROI(i).ID);
		if ~isempty(tmp)	
			for j=1:length(tmp),
				Label(tmp(j)).Nom = ROI(i).Nom_L;
			end
		end
	end
	
%end

if dd==0.0 % <- ajout le 18/07/01
	for i=size(Label):-1:2,
		Label(i).Nom = Label(i-1).Nom;
		Dist(i) = Dist(i-1);
	end
	Label(1).Nom = nn;
	Dist(1) = dd;
end
