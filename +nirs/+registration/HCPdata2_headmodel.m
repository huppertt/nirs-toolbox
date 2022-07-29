function data=HCPdata2_headmodel(SubjectsDIR)


mesh=nirs.registration.sMRI_segment(fullfile(SubjectsDIR,'mri','T1w_hires.nii.gz'));

disp('Computing wavelet model of pial surface');

% Use the wavelet methods defined in 
% Phys Med Biol. 2009 Oct 21;54(20):6383-413. doi: 10.1088/0031-9155/54/20/023. Epub 2009 Oct 7.
% Topographic localization of brain activation in diffuse optical imaging using spherical wavelets.
% Abdelnour F1, Schmidt B, Huppert TJ.
J=5;  %levels of wavelets in the model (=ico5)

options=struct('base_mesh', 'ico', 'keep_subdivision', true);
[vertex,face] = compute_semiregular_sphere(J,options);

[lsv,lsf]=read_surf(fullfile(SubjectsDIR,'surf','lh.sphere.reg'));
[rsv,rsf]=read_surf(fullfile(SubjectsDIR,'surf','rh.sphere.reg'));



kl=dsearchn(lsv/100,vertex{end}');
kr=dsearchn(rsv/100,vertex{end}');

% Create the surfaces
[lv,~]=read_surf(fullfile(SubjectsDIR,'surf','lh.pial'));
[rv,~]=read_surf(fullfile(SubjectsDIR,'surf','rh.pial'));

mesh(end+1)=nirs.core.Mesh([lv(kl,:); rv(kr,:)],[face{end}'; face{end}'+length(kl)]);  

return;