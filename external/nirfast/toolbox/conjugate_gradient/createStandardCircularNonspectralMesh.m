function [ mesh_anom mesh ] = createStandardCircularNonspectralMesh(contrast)
%CREATESTANDARDCIRCULARNONSPECTRALMESH Summary of this function goes here
%   Detailed explanation goes here
if nargin == 0
    contrast = 2;
end
mesh = load_mesh('circle2000_86_stnd');
angle = (pi/3):(2*pi/3):(2*pi);
muaBackground = mean(mesh.mua);
musBackground = mean(mesh.mus);
mua = [1 contrast contrast];
mus = [contrast contrast 1];
r_full = 43;
r_centre = r_full/2;
r_circle = 7.5;
mesh_anom = mesh;
for ii = 1:length(angle)
blob.x = r_centre*sin(angle(ii));
blob.y = r_centre*cos(angle(ii));
blob.r = r_circle;
blob.mua = mua(ii)*muaBackground;
blob.mus = mus(ii)*musBackground;
blob.region = ii;
mesh_anom = add_blob_stnd(mesh_anom,blob);
end
end

