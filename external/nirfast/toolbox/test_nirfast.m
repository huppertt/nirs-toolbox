function test_nirfast(loc)

% test_nirfast(loc)
%
% runs various test scripts for nirfast
%
% loc is the location of the nirfast folder
%   ie: test_nirfast('C:\Nirfast');

% If no input argument, detect the root of nirfast
if nargin == 0
    foo = which('nirfast');
    if isempty(foo)
        error(' Can''t find nirfast main script.')
    end
    foo = regexp(fileparts(foo),filesep,'split');
    nirfastroot = fullfile(foo{1:end-1});
    if ~ispc
        nirfastroot = [filesep nirfastroot];
    end
    loc = nirfastroot;
end
    
%% standard 2d simulation
disp('TEST_NIRFAST: standard 2d simulation');

circle2000_86_stnd = load_mesh(fullfile(loc,'meshes','standard','circle2000_86_stnd'));
blob.x=-14;
blob.y=14;
blob.r=12;
blob.z=0;
blob.mua=0.03;
blob.mus=3;
blob.region=1;
circle2000_86_stnd_anom = add_blob(circle2000_86_stnd,blob);
plotmesh(circle2000_86_stnd_anom,1);
circle2000_86_stnd_anom_data = femdata(circle2000_86_stnd_anom,100);
mesh =circle2000_86_stnd_anom;
plotimage(mesh,log(circle2000_86_stnd_anom_data.phi(:,1))); pause(0.2);
plotimage(mesh,log(circle2000_86_stnd_anom_data.phi(:,2))); pause(0.2);
plot_data(circle2000_86_stnd_anom_data);
circle2000_86_stnd_anom_data_noise = add_noise(circle2000_86_stnd_anom_data,2,2);

% diffuse
lambda.type='Automatic';
lambda.value=10;
[mesh,pj] = reconstruct(circle2000_86_stnd,[30 30],100,circle2000_86_stnd_anom_data_noise,40,lambda,'',0);
read_solution(mesh,'');

% soft priors
lambda.type='Automatic';
lambda.value=10;
[mesh,pj] = reconstruct_stnd_spatial(circle2000_86_stnd,[30 30],100,circle2000_86_stnd_anom_data_noise,40,10,'',0);
read_solution(mesh,'');

% hard priors
val.mua=0.01;
val.mus=1;
circle2000_86_stnd_anom = set_mesh(circle2000_86_stnd_anom,1,val);
lambda.type='Automatic';
lambda.value=10;
[mesh,pj] = reconstruct_stnd_region(circle2000_86_stnd_anom,100,circle2000_86_stnd_anom_data_noise,40,10,'',0,[0 1]);
read_solution(mesh,'');

disp('TEST_NIRFAST: complete');
clearvars -except loc
close all

%% fluorescence 3d simulation cw
disp('TEST_NIRFAST: fluorescence 3d simulation cw');

cylinder_fl = load_mesh(fullfile(loc,'meshes','fl','cylinder_fl'));
blob.x=0;
blob.y=-10;
blob.r=12;
blob.z=0;
blob.muaf=0.04;
blob.region=1;
cylinder_fl_anom = add_blob(cylinder_fl,blob);
cylinder_fl_anom_data = femdata(cylinder_fl_anom,0);
plot_data(cylinder_fl_anom_data);
cylinder_fl_anom_data_noise = add_noise(cylinder_fl_anom_data,2,2);
lambda.type='Automatic';
lambda.value=10;
[mesh,pj] = reconstruct(cylinder_fl,[18 18 16],0,cylinder_fl_anom_data_noise,40,lambda,'',0);
read_solution(mesh,'');

disp('TEST_NIRFAST: complete');
clearvars -except loc
close all

%% spectral 2d real
disp('TEST_NIRFAST: spectral 2d real');

[data_cal,mesh_cal]=calibrate_spectral(fullfile(loc,'data','circle_spec','homog.paa'),...
    fullfile(loc,'data','circle_spec','anom.paa'),fullfile(loc,'data','circle_spec','circle_spec'),...
    fullfile(loc,'data','circle_spec','circle_spec'),100,5);
lambda.type='Automatic';
lambda.value=10;
[mesh,pj] = reconstruct(mesh_cal,[30 30],100,data_cal,40,lambda,'',0);
read_solution(mesh,'');

disp('TEST_NIRFAST: complete');
clearvars -except loc
close all

%% create mesh simple shapes
disp('TEST_NIRFAST: create mesh simple shapes');

sizevar.xc=0;
sizevar.yc=0;
sizevar.r=33;
sizevar.dist=3;
create_mesh('','Circle',sizevar,'stnd');
clear sizevar

sizevar.xc=0;
sizevar.yc=0;
sizevar.zc=0;
sizevar.r=33;
sizevar.height=33;
sizevar.dist=3;
create_mesh('','Cylinder',sizevar,'stnd');
clear sizevar

disp('TEST_NIRFAST: complete');
clearvars -except loc
close all

%% create mesh from surface
disp('TEST_NIRFAST: create mesh from surface');

checkerboard3dmm_wrapper(fullfile(loc,'meshes','meshing examples','p1915_1.inp'),fullfile('stnd'),'stnd');

surf2nirfast_bem(fullfile(loc,'meshes','meshing examples','3026_1.inp'),fullfile('stnd_bem'));

disp('TEST_NIRFAST: complete');
clearvars -except loc
close all

%% vtk export
disp('TEST_NIRFAST: vtk export');

circle2000_86_stnd = load_mesh(fullfile(loc,'meshes','standard','circle2000_86_stnd'));
nirfast2vtk(circle2000_86_stnd,'test.vtk');

disp('TEST_NIRFAST: complete');
clearvars -except loc
close all;
