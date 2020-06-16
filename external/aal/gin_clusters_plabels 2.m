function varargout = gin_clusters_plabels(varargin)
% Display and analysis of SPM{.}
% FORMAT TabDat = gin_clusters_plabels('List',SPM,hReg,[Num,Dis,Str])
% Summary list of local maxima for entire volume of interest
% 
% SPM    - structure containing SPM, distribution & filtering details
%        - spm2
%
% (see spm_getSPM for further details of xSPM structures)
%
% hReg   - Handle of results section XYZ registry (see spm_results_ui.m)
%
% Num    - number of maxima per cluster
% Dis    - distance among clusters (mm)
% Str    - header string
%
% TabDat - Structure containing table data
%        - fields are
% .tit   - table Title (string)
% .hdr   - table header (2x3 cell array)
% .fmt   - fprintf format strings for table data (1x3 cell array)
% .str   - table filtering note (string)
% .ftr   - table footnote information (4x2 cell array)
% .dat   - table data (Nx3 cell array)
%
%                           ----------------
% See spm_list.m for more details
%
%                           ----------------
% FORMAT gin_clusters_plabels('TxtList',TabDat,c)
% Prints a tab-delimited text version of the table
% TabDat - Structure containing table data (format as above)
% c      - Column of table data to start text table at
%          (E.g. c=3 doesn't print set-level results contained in columns 1 & 2)
%
%                           ----------------
% This function loads the labelised atlas, MNI VOIs & Labels, and the distances
% between VOIs and Labels.
% For each cluster (identified by its first local maximum), 
% this function search for the intersected regions with the cluster,
% and displays the table :
%		- local maximum for the cluster
%		- the intersected regions for the cluster.
%		- the voxels pourcents for the cluster belongs to 
% 		each of the regions (labels)
% To Compute Labels & Distances we use the 'gin_det_plabels'
% routine. 
% ______________________________________________________________________
%
% gin_clusters_plabels			B Landeau 26/01/2006
% _______________________________________________________________________
% 


% satellite figure
%-----------------------------------------------------------------------
SatWindow = spm_figure('FindWin','Satellite');



%=======================================================================
switch lower(varargin{1}), case 'list'                            %-List
%=======================================================================
% FORMAT TabDat = gin_clusters_plabels('list',SPM,hReg)

%-Tolerance for p-value underflow, when computing equivalent Z's
%-----------------------------------------------------------------------
tol = eps*10;

%-Parse arguments and set maxima number and separation
%-----------------------------------------------------------------------
if nargin < 2,	error('insufficient arguments'),     end
if nargin < 3,	hReg = []; else, hReg = varargin{3}; end


%-Get current location (to highlight selected voxel in table)
%-----------------------------------------------------------------------
% xyzmm     = spm_results_ui('GetCoords');
% Modif BL 09/09/03

%-Extract data from xSPM
%-----------------------------------------------------------------------
S     = varargin{2}.S;
R     = varargin{2}.R;
FWHM  = varargin{2}.FWHM;
VOX   = varargin{2}.VOX;
n     = varargin{2}.n;
STAT  = varargin{2}.STAT;
df    = varargin{2}.df;
u     = varargin{2}.u;
M     = varargin{2}.M;
v2r   = 1/prod(FWHM(~isinf(FWHM)));			%-voxels to resels
k     = varargin{2}.k*v2r;
try
  QPs   = varargin{2}.Ps;					% Needed for FDR
  QPs   = sort(QPs(:));
catch
end


%-get number and separation for maxima to be reported
%-----------------------------------------------------------------------
if length(varargin) > 3

	Num    = varargin{4};		% number of maxima per cluster
	Dis    = varargin{5};		% distance among clusters (mm)
else
	Num    = 3;
	Dis    = 8;
end
if length(varargin) > 5
	Title  = varargin{6};
else
	Title  = ['volume summary',...
			' (labels and percentages per cluster)'];
end
% Modif BL 09/09/03


%-Setup graphics panel
%-----------------------------------------------------------------------
spm('Pointer','Watch')
if ~isempty(SatWindow)
	Fgraph = SatWindow;
	figure(Fgraph);
else
	Fgraph = spm_figure('GetWin','Graphics');
end
spm_results_ui('Clear',Fgraph)
FS    = spm('FontSizes');			%-Scaled font sizes
PF    = spm_platform('fonts');		%-Font names (for this platform)


%-Table header & footer
%=======================================================================

%-Table axes & Title
%----------------------------------------------------------------------
if ~isempty(SatWindow), ht = 0.85; bot = .14; else, ht = 0.4; bot = .1; end;
ht = 0.85; 
bot = .14;

%if STAT == 'P'
%	Title = 'Posterior Probabilities';
%end
% Modif BL 09/09/03

	
hAx   = axes('Position',[0.025 bot 0.9 ht],...
	'DefaultTextFontSize',FS(8),...
	'DefaultTextInterpreter','Tex',...
	'DefaultTextVerticalAlignment','Baseline',...
	'Units','points',...
	'Visible','off');

AxPos = get(hAx,'Position'); set(hAx,'YLim',[0,AxPos(4)])
dy    = FS(9);
y     = floor(AxPos(4)) - dy;

%text(0,y,['Statistics:  \it\fontsize{',num2str(FS(9)),'}',Title],...
%	'FontSize',FS(11),'FontWeight','Bold');	y = y - dy/2;
%line([0 1],[y y],'LineWidth',3,'Color','r'),	y = y - 9*dy/8;
% Modif BL 09/09/03

text(0,y,['Working Dir : ',varargin{2}.swd],'FontSize',FS(11),'FontWeight','Bold','Editing','on');	
y = y - 1.5*dy;
text(0,y,varargin{2}.title,'FontSize',FS(11),'FontWeight','Bold');	
y = y - 1.5*dy;
text(0,y,['Labels : \it\fontsize{',num2str(FS(9)),'}',Title],...
	'FontSize',FS(11),'FontWeight','Bold');	
y = y - dy/2;
line([0 1],[y y],'LineWidth',3,'Color','r');
y = y - 9*dy/8;


%-Construct table header
%-----------------------------------------------------------------------
set(gca,'DefaultTextFontName',PF.helvetica,'DefaultTextFontSize',FS(8))

Hc = [];

text(0.01,y - dy/2,'x,y,z {mm}','FontSize',FS(9));
text(0.15,y - dy/2,'Label','FontSize',FS(9))
text(0.40,y - dy/2,'% Cluster','FontSize',FS(9));
text(0.55,y - dy/2,'Nb Vx Cluster','FontSize',FS(9))
text(0.70,y - dy/2,'% Label','FontSize',FS(9));
text(0.85,y - dy/2,'Nb Vx Label','FontSize',FS(9));

% Modif BL 12/05/04
label='nom du label';


%-Headers for text table ... 
%-----------------------------------------------------------------------
TabDat.tit = Title;

TabDat.hdr = {	'',		'x,y,z {mm}';...
		'',		label;...
		'',		'% Cluster';...
		'',		'Nb Vx Cluster';...
		'',		'% Label';...
		'',		'Nb Vx Label'}';

TabDat.fmt = { '%3.0f %3.0f %3.0f','%s','%3.2f','%6d','%3.2f','%6d'};
% Modif BL 12/05/04


%-Column Locations
%-----------------------------------------------------------------------
%tCol       = [  0.00      0.07 ...				%-Set
%	        0.16      0.26      0.34 ...			%-Cluster
%	        0.46      0.55      0.62      0.71      0.80 ...%-Voxel
%               0.92];						%-XYZ
% Modif BL 09/09/03

% move to next vertial postion marker
%-----------------------------------------------------------------------
y     = y - 7*dy/4;
line([0 1],[y y],'LineWidth',1,'Color','r')
y     = y - 5*dy/4;
y0    = y;


%-Table filtering note
%-----------------------------------------------------------------------
if isinf(Num)
	TabDat.str = sprintf('table shows first local maximum per cluster.');
else
	%TabDat.str = sprintf(['table shows at most %s local maxima ',...
	%	'> %.1fmm apart per cluster'],Num,Dis);
	TabDat.str = sprintf('table shows first local maximum per cluster.');
end
text(0.5,4,TabDat.str,'HorizontalAlignment','Center','FontName',PF.helvetica,...
    'FontSize',FS(8),'FontAngle','Italic')
% Modif BL 09/09/03



%-Volume, resels and smoothness (if classical inference)
%-----------------------------------------------------------------------
%if STAT ~= 'P' ????????
%-----------------------------------------------------------------------
FWHMmm          = FWHM.*VOX; 				% FWHM {mm}
Pz              = spm_P(1,0,u,df,STAT,1,n,S);
Pu              = spm_P(1,0,u,df,STAT,R,n,S);
try
  Qu              = spm_P_FDR(u,df,STAT,n,QPs);
catch
  Qu = nan;
end
[P Pn Em En ] = spm_P(1,k,u,df,STAT,R,n,S);


%-Footnote with SPM parameters
%-----------------------------------------------------------------------
line([0 1],[0 0],'LineWidth',1,'Color','r')
set(gca,'DefaultTextFontName',PF.helvetica,...
	'DefaultTextInterpreter','None','DefaultTextFontSize',FS(8))
TabDat.ftr    = cell(5,2);
TabDat.ftr{1} = ...
	sprintf('Height threshold: %c = %0.2f, p = %0.3f (%0.3f)',...
		 STAT,u,Pz,Pu);
TabDat.ftr{2} = ...
	sprintf('Extent threshold: k = %0.0f voxels, p = %0.3f (%0.3f)',...
	         k/v2r,Pn,P);
TabDat.ftr{3} = ...
	sprintf('Expected voxels per cluster, <k> = %0.3f',En/v2r);
TabDat.ftr{4} = ...
	sprintf('Expected number of clusters, <c> = %0.2f',Em*Pn);
TabDat.ftr{5} = ...
	sprintf('Expected false discovery rate, <= %0.2f',Qu);
TabDat.ftr{6} = ...
	sprintf('Degrees of freedom = [%0.1f, %0.1f]',df);
TabDat.ftr{7} = ...
	sprintf(['Smoothness FWHM = %0.1f %0.1f %0.1f {mm} ',...
		 ' = %0.1f %0.1f %0.1f {voxels}'],FWHMmm,FWHM);
TabDat.ftr{8} = ...
	sprintf('Search vol: %0.0f cmm; %0.0f voxels; %0.1f resels',S*prod(VOX),S,R(end));
TabDat.ftr{9} = ...
	sprintf(['Voxel size: [%0.1f, %0.1f, %0.1f] mm ',...
		' (1 resel = %0.2f voxels)'],VOX,prod(FWHM));

text(0.0,-1*dy,TabDat.ftr{1},...
	'UserData',[u,Pz,Pu,Qu],'ButtonDownFcn','get(gcbo,''UserData'')')
text(0.0,-2*dy,TabDat.ftr{2},...
	'UserData',[k/v2r,Pn,P],'ButtonDownFcn','get(gcbo,''UserData'')')
text(0.0,-3*dy,TabDat.ftr{3},...
	'UserData',En/v2r,'ButtonDownFcn','get(gcbo,''UserData'')')
text(0.0,-4*dy,TabDat.ftr{4},...
	'UserData',Em*Pn,'ButtonDownFcn','get(gcbo,''UserData'')')
text(0.0,-5*dy,TabDat.ftr{5},...
	'UserData',Qu,'ButtonDownFcn','get(gcbo,''UserData'')')
text(0.5,-1*dy,TabDat.ftr{6},...
	'UserData',df,'ButtonDownFcn','get(gcbo,''UserData'')')
text(0.5,-2*dy,TabDat.ftr{7},...
	'UserData',FWHMmm,'ButtonDownFcn','get(gcbo,''UserData'')')
text(0.5,-3*dy,TabDat.ftr{8},...
	'UserData',[S*prod(VOX),S,R(end)],...
	'ButtonDownFcn','get(gcbo,''UserData'')')
text(0.5,-4*dy,TabDat.ftr{9},...
	'UserData',[VOX,prod(FWHM)],...
	'ButtonDownFcn','get(gcbo,''UserData'')')



%-Characterize excursion set in terms of maxima
% (sorted on Z values and grouped by regions)
%=======================================================================
if ~length(varargin{2}.Z)
	text(0.5,y-6*dy,'no suprathreshold clusters',...
		'HorizontalAlignment','Center',...
		'FontAngle','Italic','FontWeight','Bold',...
		'FontSize',FS(16),'Color',[1,1,1]*.5);
	TabDat.dat = cell(0,6);
	varargout  = {TabDat};
	spm('Pointer','Arrow')
	return
end

% Includes Darren Gitelman's code for working around
% spm_max for conjunctions with negative thresholds
%-----------------------------------------------------------------------
minz        = abs(min(min(varargin{2}.Z)));
zscores     = 1 + minz + varargin{2}.Z;
[N Z XYZ A] = spm_max(zscores,varargin{2}.XYZ);
Z           = Z - minz - 1;

%-Convert cluster sizes from voxels to resels
%-----------------------------------------------------------------------
if isfield(varargin{2},'VRvp')
	V2R = spm_get_data(varargin{2}.VRvp,XYZ);
else
	V2R = v2r;
end
N           = N.*V2R;

%-Convert maxima locations from voxels to mm
%-----------------------------------------------------------------------
XYZmm = M(1:3,:)*[XYZ; ones(1,size(XYZ,2))];

% load MNI ROIs & Labels
%-----------------------------------------------------------------------
spm_get_defaults;
global defaults;

flag_flip = 0;	% Modif BL 4/10/03
% defaults = L/R -> MNI = L/L
%	if defaults.analyze.flip 
%		defaults.analyze.flip = 0;
%		flag_flip = 1;
%	end

	%MNI=spm_get(1,'.img','select an labelised atlas');
    MNI=spm_select(1,'image','Select an labelised atlas',[],...
        spm_file(which(mfilename),'fpath'));

	MNID=spm_vol(MNI);
	MNIY = spm_read_vols(MNID);
	%eval(['load ''',strrep(MNI,'.nii','_List.mat''')]);
    load(spm_file(MNI,'suffix','_List','ext','.mat'));
	%ROI.ID
	%ROI.Nom_C
	%ROI.Nom_L

%	if flag_flip 
%		defaults.analyze.flip = 1;
%		flag_flip = 0;
%	end
		
voiCluster = M(1,1)*M(2,2)*M(3,3);
voiROI = MNID.mat(1,1)*MNID.mat(2,2)*MNID.mat(3,3);
if voiCluster ~= voiROI 
	fact = abs(voiROI/voiCluster);
else
	fact = 1;
end


%-Table proper (& note all data in cell array)
%=======================================================================

%-Pagination variables
%-----------------------------------------------------------------------
hPage = [];
set(gca,'DefaultTextFontName',PF.courier,'DefaultTextFontSize',FS(7))


%-Set-level p values {c} - do not display if reporting a single cluster
%-----------------------------------------------------------------------
c     = max(A);					%-Number of clusters
TabLin     = 1;					%-Table data line


%-Local maxima p-values & statistics
%-----------------------------------------------------------------------
HlistXYZ = [];
while sum(~isnan(Z)),

	% Paginate if necessary
	%---------------------------------------------------------------
	if y < min(Num + 1,3)*dy

		% added Fgraph term to paginate on Satellite window
		%-------------------------------------------------------
		h     = text(0.5,-5*dy,...
			sprintf('Page %d',spm_figure('#page',Fgraph)),...
			'FontName',PF.helvetica,'FontAngle','Italic',...
			'FontSize',FS(8));

		spm_figure('NewPage',[hPage,h])
		hPage = [];
		y     = y0;
	end

    	%-Find largest remaining local maximum
    	%---------------------------------------------------------------
	[U,i]   = max(Z);			% largest maxima
	j       = find(A == A(i));		% maxima in cluster

	%- MODIF BL
	%-Jump to voxel nearest current location
	[xyzmm,ii] = spm_XYZreg('NearestXYZ', XYZmm(:,i), varargin{2}.XYZmm);
	
	%save ficp xyzmm -append
	%-Find selected cluster
	B         = spm_clusters(varargin{2}.XYZ);
	jj        = find(B == B(ii));

	loc	    = varargin{2}.XYZmm(:,jj); 
	
	%-Compute Labels & Percent
    %---------------------------------------------------------------
		[Label, Perc, nbvCluster, PercROI, nbvROI] = gin_det_plabels(loc, MNIY, MNID, ROI,fact);


		h     = text(0.01,y,sprintf(TabDat.fmt{1},XYZmm(:,i)),...
			'FontWeight','Bold',...
			'Tag','ListXYZ',...
			'ButtonDownFcn',...
			'spm_mip_ui(''SetCoords'',get(gcbo,''UserData''));',...
			'Interruptible','off','BusyAction','Cancel',...
			'UserData',XYZmm(:,i));
		HlistXYZ = [HlistXYZ, h];
		hPage = [hPage, h];
 
 		for kkk=size(Perc,2):-1:1,
			%if Perc(kkk) >= 1.0
			if y < 2*dy
				h = text(0.5,-5*dy,sprintf('Page %d',...
					spm_figure('#page')),...
					'FontName',PF.helvetica,...
					'FontAngle','Italic',...
					'FontSize',FS(8));
				spm_figure('NewPage',[hPage,h])
				hPage = [];
				y     = y0;
			end

			h     = text(0.13,y,sprintf(TabDat.fmt{2},Label(kkk).Nom),	'FontWeight','Bold',...
				'UserData',Label(kkk).Nom,'ButtonDownFcn','get(gcbo,''UserData'')');
			hPage = [hPage, h];

			h     = text(0.40,y,sprintf(TabDat.fmt{3},Perc(kkk)),	'FontWeight','Bold',...
				'UserData',Perc(kkk),'ButtonDownFcn','get(gcbo,''UserData'')');
			hPage = [hPage, h];
			
			h     = text(0.55,y,sprintf(TabDat.fmt{4},nbvCluster(kkk)),	'FontWeight','Bold',...
				'UserData',nbvCluster(kkk),'ButtonDownFcn','get(gcbo,''UserData'')');
			hPage = [hPage, h];

			h     = text(0.70,y,sprintf(TabDat.fmt{5},PercROI(kkk)),	'FontWeight','Bold',...
				'UserData',PercROI(kkk),'ButtonDownFcn','get(gcbo,''UserData'')');
			hPage = [hPage, h];

			h     = text(0.85,y,sprintf(TabDat.fmt{6},nbvROI(kkk)),	'FontWeight','Bold',...
				'UserData',nbvROI(kkk),'ButtonDownFcn','get(gcbo,''UserData'')');
			hPage = [hPage, h];

			y = y - dy;
			[TabDat.dat{TabLin,1:6}] = deal(XYZmm(:,i),Label(kkk).Nom,Perc(kkk),...
				nbvCluster(kkk),PercROI(kkk),nbvROI(kkk));
			TabLin = TabLin + 1;

			%end
		end
    %---------------------------------------------------------------

	
	Z(j) = NaN;		% Set local maxima to NaN
end				% end region


%-Number and register last page (if paginated)
%-Changed to use Fgraph for numbering
%-----------------------------------------------------------------------
if spm_figure('#page',Fgraph)>1
	h = text(0.5,-5*dy,sprintf('Page %d/%d',spm_figure('#page',Fgraph)*[1,1]),...
		'FontName',PF.helvetica,'FontSize',FS(8),'FontAngle','Italic');
	spm_figure('NewPage',[hPage,h])
end

%-End: Store TabDat in UserData of axes & reset pointer
%=======================================================================
h      = uicontextmenu('Tag','TabDat',...
		'UserData',TabDat);
set(gca,'UIContextMenu',h,...
	'Visible','on',...
	'XColor','w','YColor','w')
uimenu(h,'Label','Table')
uimenu(h,'Separator','on','Label','Print text table',...
	'Tag','TD_TxtTab',...
	'CallBack',...
	'gin_clusters_plabels(''txtlist'',get(get(gcbo,''Parent''),''UserData''),3)',...
	'Interruptible','off','BusyAction','Cancel');
uimenu(h,'Separator','off','Label','Extract table data structure',...
	'Tag','TD_Xdat',...
	'CallBack','get(get(gcbo,''Parent''),''UserData'')',...
	'Interruptible','off','BusyAction','Cancel');
uimenu(h,'Separator','on','Label','help',...
	'Tag','TD_Xdat',...
	'CallBack','spm_help(''gin_clusters_plabels'')',...
	'Interruptible','off','BusyAction','Cancel');

%-Setup registry
%-----------------------------------------------------------------------
set(hAx,'UserData',struct('hReg',hReg,'HlistXYZ',HlistXYZ))
%spm_XYZreg('Add2Reg',hReg,hAx,'spm_list');

%-Return TabDat structure & reset pointer
%-----------------------------------------------------------------------
varargout = {TabDat};
spm('Pointer','Arrow')



%=======================================================================
case 'txtlist',                                  %-Print ASCII text table
%=======================================================================
% FORMAT gin_clusters_plabels('TxtList',TabDat,c)

if nargin<2, error('Insufficient arguments'), end
if nargin<3, c=1; else, c=varargin{3}; end
TabDat = varargin{2};

%-Table Title
%-----------------------------------------------------------------------
fprintf('\n\nSTATISTICS: %s\n',TabDat.tit)
fprintf('%c','='*ones(1,80)), fprintf('\n')

%-Table header
%-----------------------------------------------------------------------
fprintf('%s\t',TabDat.hdr{1,c:end-1}), fprintf('%s\n',TabDat.hdr{1,end})
fprintf('%s\t',TabDat.hdr{2,c:end-1}), fprintf('%s\n',TabDat.hdr{2,end})
fprintf('%c','-'*ones(1,80)), fprintf('\n')

%-Table data
%-----------------------------------------------------------------------
for i = 1:size(TabDat.dat,1)
	for j=c:size(TabDat.dat,2)
		fprintf(TabDat.fmt{j},TabDat.dat{i,j})
		fprintf('\t')
	end
	fprintf('\n')
end
for i=1:max(1,11-size(TabDat.dat,1)), fprintf('\n'), end
fprintf('%s\n',TabDat.str)
fprintf('%c','-'*ones(1,80)), fprintf('\n')

%-Table footer
%-----------------------------------------------------------------------
fprintf('%s\n',TabDat.ftr{:})
fprintf('%c','='*ones(1,80)), fprintf('\n\n')



%=======================================================================
otherwise                                        %-Unknown action string
%=======================================================================
error('Unknown action string')


%=======================================================================
end
%=======================================================================
