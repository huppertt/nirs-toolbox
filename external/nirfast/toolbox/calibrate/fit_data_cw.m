function [mua,mus,lnI_offset,fem_data] = fit_data_cw(mesh,data,iteration,nographs)
                                                           
% [mua,mus,lnI_offset,fem_data] = fit_data(mesh,...
%                                   data,iteration)
%
% fits data to a given model to give initial guesses for
% reconstruction as well as data offsets
%
% mesh is the input mesh (variable)
% data is the boundary data (variable)
% iteration is the number of iterations for fitting
% outputs are the intial guesses
% nographs is a flag for displaying graphs


if ~exist('nographs','var')
    nographs = 0;
end

frequency = 0;
if size(data,2) > 1
    data = data(:,1);
end

% calculate the source / detector distance for each combination.
% also make sure that proper sources / detectors are used
dist = zeros(length(mesh.link),1);
for i = 1:length(mesh.link)
    snum = mesh.link(i,1);
    mnum = mesh.link(i,2);
    snum = mesh.source.num == snum;
    mnum = mesh.meas.num == mnum;
    if sum(snum)==0 || sum(mnum)==0
        dist(i,1)=0;
        mesh.link(i,3)=0;
    else
        dist(i,1) = sqrt(sum((mesh.source.coord(snum,:) - ...
        mesh.meas.coord(mnum,:)).^2,2)); 
    end
end

% get an index from link file of data to actually use
linki = logical(mesh.link(:,3));

% Set lnrI, lnr and phase!
lnrI = log(data(:,1).*dist);
lnI = log(data(:,1));

if nographs == 0
    figure;
    subplot(2,2,1)
    plot(dist(linki),lnrI(linki),'.')
    ylabel('lnrI');
    xlabel('Source / Detector distance');
    drawnow
    pause(0.001)
end

% Calculate the coeff of a polynomial fit of distance vs. lnrI
m1 = polyfit(dist(linki),lnrI(linki),1); m1 = m1(1);

% fit data using an analytical model
% based on Pogue paper
omega = 2*pi*frequency*1e6;
c=(3e11./mean(mesh.ri));

mua = 0.01; mus = mean(mesh.mus); kappa = 1/(3*(mua+mus));

for i = 1 : 25
  f1 = ((mua^2 + (omega/c)^2) / (kappa^2))^0.25;
  alpha = -f1*cos(0.5*atan2(omega/c,mua));
  
  mua = mua + 0.0001; 
  kappa = 1/(3*(mua+mus));
  f1 = ((mua^2 + (omega/c)^2) / (kappa^2))^0.25;
  alpha1 = -f1*cos(0.5*atan2(omega/c,mua));
  mua = mua - 0.0001;   
  
  da = (alpha-alpha1)/0.0001;
  mua = mua - (m1-alpha)/da;
  mua = abs(mua);
  kappa = 1/(3*(mua+mus));
end

disp('Global values calculated from Analytical fit');
disp(['Absorption = ' num2str(mua) ' mm-1']);
disp('============================================');

% Set the global values onto mesh.
mesh.mua(:) = mua;
mesh.kappa = 1./(3*(mesh.mua+mesh.mus));

% Fit for mua and mus using FEM
% dist = dist_orig;
jj = 0;
while jj ~= iteration
  [fem_data]=femdata(mesh,frequency);
  fem_data = fem_data.paa;
    
  femlnrI = log(fem_data(:,1).*dist);
  
  alpha0 = polyfit(dist(linki),femlnrI(linki),1); alpha0 = alpha0(1);
  
  mesh.mua(:) = mesh.mua(:)+0.0001;
  mesh.kappa = 1./(3*(mesh.mua+mesh.mus));
  
  [fem_data]=femdata(mesh,frequency);
  fem_data = fem_data.paa;
  
  femlnrI = log(fem_data(:,1).*dist);
  
  alpha1 = polyfit(dist(linki),femlnrI(linki),1); alpha1 = alpha1(1);
  
  mesh.mua(:) = mesh.mua(:)-0.0001;
  mesh.kappa = 1./(3*(mesh.mua+mesh.mus));
  
  da = (alpha0-alpha1)/0.0001;
  
  mesh.mua(:) = mesh.mua(:) - (m1-alpha0)/da;
  mesh.kappa = 1./(3*(mesh.mua+mesh.mus));
  
  err_a = abs(mean(mesh.mua)-mua)./mua;
  
  if ((err_a < 0.001))
    jj = iteration;
  else
    jj = jj + 1;
  end
  
  mua = mean(mesh.mua);
  mus = mean(mesh.mus);
  if mua < 0 || mus < 0
    errordlg('Negative mua or mus was calculated. This may be caused by bad/noisy data.','NIRFAST Error');
    error('Negative mua or mus was calculated. This may be caused by bad/noisy data.');
  end
  
end

disp('Global values calculated from Numerical fit');
  disp(['Absorption = ' num2str(mua) ' mm-1 with error of ' num2str(err_a)]);
  disp('-------------------------------------------------');


% Calculate data based on these global values
[fem_data]=femdata(mesh,frequency);
fem_data = fem_data.paa;


% Arrange data to calculate offset
femlnI = log(fem_data(:,1));

% Find offset
% make sure any original NaNs are present
lnI = log(data(:,1));

% Set offset based on particular source / detector
% we do this because data is not symmetrical!
n = max(mesh.link(:,1));
lnI_offset = zeros(size(lnI,1),1);
for i=1:n
    source_all = mesh.link(:,1)==i;
    source_used = and(mesh.link(:,1)==i , linki);
    lnI_offset(source_all) = mean(lnI(source_used)-femlnI(source_used));
end

if nographs == 0
    subplot(2,2,3:4);
    plot(lnI(linki)-lnI_offset(linki),'k');
    hold on
    plot(femlnI(linki),'r--');
    axis tight;
    xlabel('log Amplitude');
    legend('Measured-Offset','Model');
end

end
