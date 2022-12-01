function mesh=HCPdata2_headmodel(subjdir)

T1file=rdir(fullfile(subjdir,'T1w','*','mri','T1w_hires.nii.gz'));
mesh=nirs.registration.sMRI_segment(T1file(1).name);

disp('Computing wavelet model of pial surface');

% Use the wavelet methods defined in 
% Phys Med Biol. 2009 Oct 21;54(20):6383-413. doi: 10.1088/0031-9155/54/20/023. Epub 2009 Oct 7.
% Topographic localization of brain activation in diffuse optical imaging using spherical wavelets.
% Abdelnour F1, Schmidt B, Huppert TJ.
J=5;  %levels of wavelets in the model (=ico5)

options=struct('base_mesh', 'ico', 'keep_subdivision', true);
[vertex,face] = compute_semiregular_sphere(J,options);

SubjectsDIR=rdir(fullfile(subjdir,'T1w','*','surf'));
[lsv,~]=read_surf(fullfile(SubjectsDIR(1).folder,'lh.sphere.reg'));
[rsv,~]=read_surf(fullfile(SubjectsDIR(1).folder,'rh.sphere.reg'));
kl=dsearchn(lsv/100,vertex{end}');
kr=dsearchn(rsv/100,vertex{end}');


gr=rdir(fullfile(subjdir,'MNINonLinear','fsaverage_LR32k','*.R.sphere.32k_fs_LR.surf.gii'));
gl=rdir(fullfile(subjdir,'MNINonLinear','fsaverage_LR32k','*.L.sphere.32k_fs_LR.surf.gii'));
gr=gifti(gr(end).name);
gl=gifti(gl(end).name);


kl2=dsearchn(gl.vertices/100,lsv(kl,:));
kr2=dsearchn(gr.vertices/100,rsv(kr,:));

% Create the surfaces
pr=rdir(fullfile(subjdir,'MNINonLinear','fsaverage_LR32k','*.R.pial_MSMAll.32k_fs_LR.surf.gii'));
pl=rdir(fullfile(subjdir,'MNINonLinear','fsaverage_LR32k','*.L.pial_MSMAll.32k_fs_LR.surf.gii'));
pr=gifti(pr(end).name);
pl=gifti(pl(end).name);

lv=pl.vertices;
rv=pr.vertices;

%[lv,~]=read_surf(fullfile(SubjectsDIR(1).folder,'lh.pial'));
%[rv,~]=read_surf(fullfile(SubjectsDIR(1).folder,'rh.pial'));


taskfile=rdir(fullfile(subjdir,'MNINonLinear','Results','*',['*_hp200_s*level2_MSMAll.feat'],['*_MSMAll.dscalar.nii']))

c=ft_read_cifti(taskfile(end).name);
lstl=find(c.brainstructure==1);
lstr=find(c.brainstructure==2);

mesh(end+1)=nirs.core.Mesh([pl.vertices(lstl(kl2),:); pr.vertices(lstr(kr2)-max(lstl),:)],[face{end}'; face{end}'+length(vertex{end})]);


T =[ 0.9964    0.0178    0.0173   -0.0000
   -0.0169    0.9957   -0.0444   -0.0000
   -0.0151    0.0429    1.0215    0.0000
   -0.4232  -17.5022   11.6967    1.0000];
mesh(end)=nirs.registration.applytform(mesh(end),inv(T));

return;