function [mua, mus, kappa, E] = calc_mua_mus(mesh,wv_array)

% [mua, mus, kappa, E] = calc_mua_mus(mesh,wv_array)
%
% Given a spectral mesh, and specific wavelengths, 
% calculates mua, mus and kappa at these wavelengths


wv_array = sort(wv_array);

%****************************************************************
% calculate absorption coefficients
E = []; mus = [];

% for mus, wavelength must be in micrometers.
wv_act = wv_array/1000;

for i = 1:length(wv_array)
  index = find(mesh.wv == wv_array(i));
  E = [E;mesh.excoef(index,:)];
  mus = [mus mesh.sa.*wv_act(i).^(-mesh.sp)];
end

mua = E*mesh.conc';
mua = mua';
kappa = 1./(3*(mus+mua));
