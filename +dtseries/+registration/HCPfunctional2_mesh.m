function data=HCPfunctional2_mesh(subjdir,task)

smIdx=4; % smoothing index

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
[lv,~]=read_surf(fullfile(SubjectsDIR(1).folder,'lh.pial'));
[rv,~]=read_surf(fullfile(SubjectsDIR(1).folder,'rh.pial'));


taskfile=rdir(fullfile(subjdir,'MNINonLinear','Results',task,[task '_hp200_s' num2str(smIdx) '_level2.feat'],['*_' task '_level2_hp200_s' num2str(smIdx) '.dscalar.nii']))

c=ft_read_cifti(taskfile(end).name);
lstl=find(c.brainstructure==1);
lstr=find(c.brainstructure==2);
flds=fields(c);
data=struct;
for i=1:length(flds)
    if(length(c.(flds{i}))==length(c.brainstructure))
        data=setfield(data,flds{i},c.(flds{i})([lstl(kl2); lstr(kr2)]));
    end
end