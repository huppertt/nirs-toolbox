function [mua,mus,lnI_offset,phase_offset,fem_data] = fit_data(mesh,...
                                                               data,...
                                                               frequency,...
                                                               iteration,...
                                                               nographs)
                                                           
% [mua,mus,lnI_offset,phase_offset,fem_data] = fit_data(mesh,...
%                                   data,frequency,iteration,nographs)
%
% fits data to a given model to give initial guesses for
% reconstruction as well as data offsets
%
% mesh is the input mesh (variable)
% data is the boundary data (variable)
% frequency is the modulation frequency (MHz)
% iteration is the number of iterations for fitting
% outputs are the intial guesses
% nographs is a flag for if the graphs are displayed


if ~exist('nographs','var')
    nographs = 0;
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
phase = data(:,2);

if nographs == 0
    figure;
    subplot(2,2,1);
    plot(dist(linki),lnrI(linki),'.')
    ylabel('lnrI');
    xlabel('Source / Detector distance');
    subplot(2,2,2);
    plot(dist(linki),phase(linki),'.')
    ylabel('Phase');
    xlabel('Source / Detector distance');
    drawnow
    pause(0.001)
end

% Calculate the coeff of a polynomial fit of distance vs. Phase or lnrI
m0 = polyfit(dist(linki),phase(linki),1); m0 = m0(1);
m1 = polyfit(dist(linki),lnrI(linki),1); m1 = m1(1);

% fit data using an analytical model
% based on Pogue paper
omega = 2*pi*frequency*1e6;
c=(3e11./mean(mesh.ri));

mua = 0.01; mus = 1; kappa = 1/(3*(mua+mus));

for i = 1 : 25
  f1 = ((mua^2 + (omega/c)^2) / (kappa^2))^0.25;
  alpha = -f1*cos(0.5*atan2(omega/c,mua));
  phi = f1*sin(0.5*atan2(omega/c,mua))*180/pi;
  
  mua = mua + 0.0001; 
  kappa = 1/(3*(mua+mus));
  f1 = ((mua^2 + (omega/c)^2) / (kappa^2))^0.25;
  alpha1 = -f1*cos(0.5*atan2(omega/c,mua));
  phi1 = f1*sin(0.5*atan2(omega/c,mua))*180/pi;
  mua = mua - 0.0001;   
  da = (alpha-alpha1)/0.0001;
  ds = (phi-phi1)/0.0001;
  mua = mua - (m1-alpha)/da;
  mus = mus - (m0-phi)/ds*0.1;
  kappa = 1/(3*(mua+mus));
  f1 = ((mua^2 + (omega/c)^2) / (kappa^2))^0.25;
  alpha = -f1*cos(0.5*atan2(omega/c,mua));
  phi = f1*sin(0.5*atan2(omega/c,mua))*180/pi;
  
  mus = mus + 0.001; 
  kappa = 1/(3*(mua+mus));
  f1 = ((mua^2 + (omega/c)^2) / (kappa^2))^0.25;
  alpha2 = -f1*cos(0.5*atan2(omega/c,mua));
  phi2 = f1*sin(0.5*atan2(omega/c,mua))*180/pi;
  mus = mus - 0.001; 
  da = (alpha-alpha2)/0.001;
  ds = (phi-phi2)/0.001;
  mua = mua - (m1-alpha)/da*0.01;
  mus = mus - (m0-phi)/ds*0.2;
  kappa = 1/(3*(mua+mus));
end

disp('Global values calculated from Analytical fit');
disp(['Absorption = ' num2str(mua) ' mm-1']);
disp(['Scatter    = ' num2str(mus) ' mm-1']);
disp('============================================');

% Set the global values onto mesh.
mesh.mua(:) = mua;
mesh.mus(:) = mus;
mesh.kappa(:) = kappa;

% Fit for mua and mus using FEM
% dist = dist_orig;
jj = 0;
while jj ~= iteration
  [fem_data]=femdata(mesh,frequency);
  fem_data = fem_data.paa;
    
  femlnrI = log(fem_data(:,1).*dist);
  femphase = fem_data(:,2);
  
  phi0 = polyfit(dist(linki),femphase(linki),1); phi0 = phi0(1);
  alpha0 = polyfit(dist(linki),femlnrI(linki),1); alpha0 = alpha0(1);
  
  mesh.mua(:) = mesh.mua(:)+0.0001;
  mesh.kappa = 1./(3*(mesh.mua+mesh.mus));
  
  [fem_data]=femdata(mesh,frequency);
  fem_data = fem_data.paa;
  
  femlnrI = log(fem_data(:,1).*dist);
  femphase = fem_data(:,2);
  
  phi1 = polyfit(dist(linki),femphase(linki),1); phi1 = phi1(1);
  alpha1 = polyfit(dist(linki),femlnrI(linki),1); alpha1 = alpha1(1);
  
  mesh.mua(:) = mesh.mua(:)-0.0001;
  mesh.kappa = 1./(3*(mesh.mua+mesh.mus));
  
  da = (alpha0-alpha1)/0.0001;
  ds = (phi0-phi1)/0.0001;
  
  mesh.mua(:) = mesh.mua(:) - (m1-alpha0)/da;
  mesh.mus(:) = mesh.mus(:) - (m0-phi0)/ds.*0.1;
  mesh.kappa = 1./(3*(mesh.mua+mesh.mus));
 [mesh.mua(1) mesh.mus(1)];
 
 %
 
 [fem_data]=femdata(mesh,frequency);
  fem_data = fem_data.paa;
  
  femlnrI = log(fem_data(:,1).*dist);
  femphase = fem_data(:,2);

  phi0 = polyfit(dist(linki),femphase(linki),1); phi0 = phi0(1);
  alpha0 = polyfit(dist(linki),femlnrI(linki),1); alpha0 = alpha0(1);
  
  mesh.mus(:) = mesh.mus(:)+0.001;
  mesh.kappa = 1./(3*(mesh.mua+mesh.mus));
  
  [fem_data]=femdata(mesh,frequency);
  fem_data = fem_data.paa;
  
  femlnrI = log(fem_data(:,1).*dist);
  femphase = fem_data(:,2);
  
  phi2 = polyfit(dist(linki),femphase(linki),1); phi2 = phi2(1);
  alpha2 = polyfit(dist(linki),femlnrI(linki),1); alpha2 = alpha2(1);
  
  mesh.mus(:) = mesh.mus(:)-0.001;
  mesh.kappa = 1./(3*(mesh.mua+mesh.mus));
  
  da = (alpha0-alpha2)/0.001;
  ds = (phi0-phi2)/0.001;
  
  mesh.mua(:) = mesh.mua(:) - abs((m1-alpha0)/da*0.002);
  mesh.mua(find(mesh.mua(:)<0)) = mesh.mua(find(mesh.mua(:)<0)) + da;
  mesh.mus(:) = mesh.mus(:) - (m0-phi0)/ds.*0.2;
  mesh.kappa = 1./(3*(mesh.mua+mesh.mus));
  
  err_a = abs(mean(mesh.mua)-mua)./mua;
  err_s = abs(mean(mesh.mus)-mus)./mus;
  
  if ((err_a < 0.001) & (err_s < 0.001))
    jj = iteration;
  else
    jj = jj + 1;
  end
  
  mua = mean(mesh.mua);
  mus = mean(mesh.mus);
  
end

disp('Global values calculated from Numerical fit');
  disp(['Absorption = ' num2str(mua) ' mm-1 with error of ' num2str(err_a)]);
  disp(['Scatter    = ' num2str(mus) ' mm-1 with error of ' num2str(err_s)]);
  disp('-------------------------------------------------');

% Calculate data based on these global values
[fem_data]=femdata(mesh,frequency);
fem_data = fem_data.paa;


% Arrange data to calculate offset
femlnI = log(fem_data(:,1));
femphase = fem_data(:,2);

% Find offset
% make sure any original NaNs are present
lnI = log(data(:,1));
phase = data(:,2);

% Set offset based on particular source / detector
% we do this because data is not symmetrical!
n = max(mesh.link(:,1));
lnI_offset = zeros(size(lnI,1),1);
phase_offset = zeros(size(phase,1),1);
for i=1:n
    source_all = mesh.link(:,1)==i;
    source_used = and(mesh.link(:,1)==i , linki);
    lnI_offset(source_all) = mean(lnI(source_used)-femlnI(source_used));
    phase_offset(source_all) = mean(phase(source_used)-femphase(source_used));
end

if nographs == 0
    subplot(2,2,3);
    plot(lnI(linki)-lnI_offset(linki),'k');
    hold on
    plot(femlnI(linki),'r--');
    axis tight;
    xlabel('log Amplitude');
    legend('Measured-Offset','Model');
    subplot(2,2,4);
    plot(phase(linki)-phase_offset(linki),'k');
    hold on
    plot(femphase(linki),'r--');
    axis tight;
    xlabel('Phase');
    legend('Measured-Offset','Model');
end

end
