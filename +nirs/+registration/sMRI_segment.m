function mesh = sMRI_segment(mfile)
%mfile='HCP201_3T_T1w_MPR1.nii.gz';

global defaults;
spm_defaults;

if(isa(mfile,'char'))
    mri=ft_read_mri(mfile);
else
    mri=mfile;
end


cfg.spmversion='spm12';
cfg.output = {'brain' 'scalp' 'skull'};


cfg.scalpsmooth    = 'no';

%     
%     cfg.method='spm';
%     cfg.nonlinear='no';
%     cfg.keepintermediate='yes';
%     [mri2] = ft_volumenormalise(cfg, mri);
%     
   
    
    

mri.coordsys='spm';


cfg.output = {'tpm' 'all'};
[segmented] = ft_volumesegment(cfg, mri);

cfg.output = {'scalp','skull','brain'};
[segmented2] = ft_volumesegment(cfg,segmented);

seg = zeros(size(segmented2.skull));
seg(segmented2.scalp)=1;
seg(segmented2.skull)=2;
seg(segmented2.brain)=3;
% 
% lst=find(segmented2.brain  & segmented.csf>max(segmented.gray,segmented.white));
% seg(lst)=4;
% 
% lst=find(segmented2.brain  & segmented.gray>max(segmented.csf,segmented.white));
% seg(lst)=5;
% 
% lst=find(segmented2.brain  & segmented.white>max(segmented.gray,segmented.csf));
% seg(lst)=3;

mesh(1)=nirs.core.Mesh;
mesh(1).transparency=.1;
mesh(2)=nirs.core.Mesh;
mesh(2).transparency=.1;
mesh(3)=nirs.core.Mesh;
mesh(3).transparency=1;


[mesh(1).faces,mesh(1).nodes]=isosurface(seg>0);
[mesh(2).faces,mesh(2).nodes]=isosurface(seg>1);
[mesh(3).faces,mesh(3).nodes]=isosurface(seg>2);


for i=1:3
    n=mesh(i).nodes;
    n(:,4)=1;
    n=n*(mri.hdr.tkrvox2ras)';
    mesh(i).nodes=n(:,1:3);
end

tbl=nirs.util.list_1020pts('?');
Pos =[-tbl.Z tbl.X tbl.Y];

[TR, TT] = icp(mesh(1).nodes',Pos');
Pos=(TR*Pos'+TT*ones(1,size(Pos,1)))';

k=dsearchn(mesh(1).nodes,Pos);
p=mesh(1).nodes(k,:);
p(:,4)=1;
Pos(:,4)=1;
Pos=Pos*(Pos\p);

fidtbl=table(tbl.Name,Pos(:,1),Pos(:,2),Pos(:,3),repmat({'10-20'},length(tbl.Name),1),...
    repmat({'mm'},length(tbl.Name),1),repmat(true,length(tbl.Name),1),...
    'VariableNames',mesh(1).fiducials.Properties.VariableNames);

mesh(1).fiducials=fidtbl;

end

