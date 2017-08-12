function [fwd_mesh,pj_error] = reconstruct_spectral_bem(fwd_mesh,...
                                                    frequency,...
                                                    data_fn,...
                                                    iteration,...
                                                    lambda,...
                                                    output_fn,...
                                                    filter_n,...
                                                    wv_array)

% [fwd_mesh,pj_error] = reconstruct_spectral_bem(fwd_mesh,...
%                                            frequency,...
%                                            data_fn,...
%                                            iteration,...
%                                            lambda,...
%                                            output_fn,...
%                                            filter_n,...
%                                            wv_array)
%                                                
% Spectrally constrained or "direct spectral" reconstruction program - bem
%
% fwd_mesh is the input mesh (variable or filename)
% frequency is the modulation frequency (MHz)
% data_fn is the data filename
% iteration is the max number of iterations
% lambda is the initial regularization value
% output_fn is the root output filename
% filter_n is the number of mean filters
% wv_array is optional wavelength array


tmainloop = tic;

% error checking
if frequency < 0
    errordlg('Frequency must be nonnegative','NIRFAST Error');
    error('Frequency must be nonnegative');
end

% initialize projection error
pj_error = [];

% If not a workspace variable, load mesh
if ischar(fwd_mesh)== 1
  fwd_mesh = load_mesh(fwd_mesh);
end
if ~strcmp(fwd_mesh.type,'spec_bem')
    errordlg('Mesh type is incorrect','NIRFAST Error');
    error('Mesh type is incorrect');
end
if exist('wv_array') == 0
  wv_array = fwd_mesh.wv;
end

%*******************************************************
% read data - This is the calibrated experimental data or simulated data
disp('Loading data and wavelength information')
data = load_data(data_fn,wv_array);
% if specified wavelength not available, terminate.
if isempty(data) || ~isfield(data,'paa')
    errordlg('Data not found or not properly formatted','NIRFAST Error');
    error('Data not found or not properly formatted');
end
[n,m] = size(data.link);  [nd,md] = size(data.paa);
data_link = data.link;

if 2*(m-2) ~= md
    errordlg('data.link does not equal data.paa','NIRFAST Error');
    error('data.link does not equal data.paa');
end

% we need log amplitude and phase in radians
anom_a = [];
anom_p = [];
k=1;
for i = 1:2:md
    data.paa(:,i) = log(data.paa(:,i));
    foo = data.paa(:,i+1)/180.0*pi;
    foo(foo<0) = foo(foo<0) + (2*pi);
    foo(foo>(2*pi)) = foo(foo>(2*pi)) - (2*pi);
    data.paa(:,i+1) = foo; clear foo
    linki = logical(data_link(:,k+2));
    anom_a = [anom_a; data.paa(linki,i)];
    anom_p = [anom_p; data.paa(linki,i+1)];
    k = k+1;
end
anom = zeros(length(anom_a)*2,1);
anom(1:2:end) = anom_a;
anom(2:2:end) = anom_p;
clear data anom_a anom_p

% check to ensure wv_array wavelengths match the wavelength list fwd_mesh
for i = 1:length(wv_array)
  tmp = find(fwd_mesh.wv == wv_array(i));
  if isempty(tmp)
      flag(i) = 0;
  else
      flag(i) = tmp(1);
  end
end
tmp = find(flag==0);
if isempty(tmp) ~= 1
    for i = 1 : length(tmp)
        disp(['ERROR: wv_array contains ' num2str(wv_array(tmp(i))) ...
              'nm which is not present in ' fwd_mesh.name,'.excoef']);
    end
    return
end
clear tmp flag i

nwv = length(wv_array);

% extinction coeff for chosen wavelengths
[junk1,junk2,junk3,E] = calc_mua_mus(fwd_mesh,wv_array);
clear junk*

%************************************************
% Initiate log file
fid_log = fopen([output_fn,'.log'],'w');
cl = fix(clock);
write_log(fid_log,'Started on %s ',date,1);
write_log(fid_log,'at %d:%d:%d\n',cl(4:6),1);
fprintf(fid_log,'Forward Mesh       = %s\n',fwd_mesh.name);
fprintf(fid_log,'Frequency          = %f MHz\n',frequency);
if ischar(data_fn) ~= 0
    fprintf(fid_log,'Data file          = %s\n',data_fn);
end
if isstruct(lambda)
    fprintf(fid_log,'Initial Regularization  = %d\n',lambda.value);
else
    fprintf(fid_log,'Initial Regularization  = %d\n',lambda);
end
fprintf(fid_log,'Filter             = %d\n',filter_n);
fprintf(fid_log,'Wavelengths Used   = ');
for i = 1 : nwv
  fprintf(fid_log,'%s ',num2str(fwd_mesh.wv(i)));
end
fprintf(fid_log,'\n\n');
fprintf(fid_log,'********************************\n');
fprintf(fid_log,'Output files:\n');

% Get string names for Chromophore files
all_sol = char(fwd_mesh.chromscattlist);
[n_allsol,junk] = size(all_sol);
str_label = [];
for i = 1:n_allsol
  tmp_string = [strcat(output_fn,'_',all_sol(i,:),'.sol')];
  str_label = [str_label tmp_string];
  fprintf(fid_log,'%s\n',tmp_string);
  clear tmp_string
end

% log initial guesses
for i = 1:n_allsol
    if strcmp(all_sol(i,1:3),'S-A')
        fprintf(fid_log,['Initial Guess ' all_sol(i,:) ' = %d\n'],...
            fwd_mesh.sa(1));
    elseif strcmp(all_sol(i,1:3),'S-P')
        fprintf(fid_log,['Initial Guess ' all_sol(i,:) ' = %d\n'],...
            fwd_mesh.sp(1));
    else
        fprintf(fid_log,['Initial Guess ' all_sol(i,:) ' = %d\n'],...
            fwd_mesh.conc(1,i));
    end
end

% start non-linear itertaion image reconstruction part
for it = 1:iteration
  
  % Compute Jacobian and ref data from each wavelength data
  disp('---------------------------------');
  disp('Building Jacobian using jacobian_spectral')
  [J,data,fwd_mesh] = jacobian_spectral_bem(fwd_mesh,frequency,wv_array);
  
  % we need log amplitude and phase in radians    
    ref_a = [];
    ref_p = [];
    k=1;
    for i = 1:2:md
        data.paa(:,i) = log(data.paa(:,i));
        foo = data.paa(:,i+1)/180.0*pi;
        foo(foo<0) = foo(foo<0) + (2*pi);
        foo(foo>(2*pi)) = foo(foo>(2*pi)) - (2*pi);
        data.paa(:,i+1) = foo; clear foo
        linki = logical(data_link(:,k+2));
        ref_a = [ref_a; data.paa(linki,i)];
        ref_p = [ref_p; data.paa(linki,i+1)];
        k = k+1;
    end
    ref = zeros(length(ref_a)*2,1);
    ref(1:2:end) = ref_a;
    ref(2:2:end) = ref_p;
    clear data ref_a ref_p
  
  % Calculate data difference:
  data_diff = (anom - ref);
  
  %********************************************************************* 
  % Update pj_error and check stopping criteria
  pj_error = [pj_error sum(abs(data_diff.^2))];  
  
  disp(['Iteration Number          = ' num2str(it)]);
  disp(['Projection error          = ' num2str(pj_error(end))]);
  
  fprintf(fid_log,'---------------------------------\n');
  fprintf(fid_log,'Iteration Number          = %d\n',it);
  fprintf(fid_log,'Projection error          = %f\n',pj_error(end));
  
  if it ~= 1
    p = (pj_error(end-1)-pj_error(end))*100/pj_error(end-1);
    disp(['Projection error change   = ' num2str(p) '%']);
    fprintf(fid_log,'Projection error change   = %f %%\n',p);
    if (p) <= 0.5 || (pj_error(it) < (10^-18))
      disp('---------------------------------');
      disp('STOPPING CRITERIA REACHED');
      fprintf(fid_log,'---------------------------------\n');
      fprintf(fid_log,'STOPPING CRITERIA REACHED\n');
      break
    end
  end
  
  %*********************************************************************
  % reduce or increase lambda based on last pj_error
  if it ~= 1
    if ((sum(data_diff.^2) < pj_error(end-1)))
      lambda = lambda./10^0.25;
    else
      lambda = lambda.*10^0.125;
    end
  end
  
  %*********************************************************************
  % Normalize Jacobian w.r.t different SI units
  conc_temp = reshape(fwd_mesh.conc, numel(fwd_mesh.conc), 1);
  N = [mean(fwd_mesh.conc) ...
       mean(fwd_mesh.sa) ...
       mean(fwd_mesh.sp)];
  nr = length(fwd_mesh.sp);
  for i = 1 : length(N)
    J(:,(i-1)*nr+1:i*nr) = J(:,(i-1)*nr+1:i*nr).*N(i);
  end
  
  % build hessian
  [nrow,ncol]=size(J);
  
  Hess = zeros(ncol);
  disp('Calculating Hessian');
  Hess = (J'*J);

  % Add regularisation
  nr = length(fwd_mesh.sp);
  nc = length(Hess)/nr;
  reg = [];
  reg_temp = eye(nr);
  for i = 1 : nc
    reg(i) = lambda.*max(diag(Hess((i-1)*nr+1:i*nr,(i-1)*nr+1:i*nr)));
    Hess((i-1)*nr+1:i*nr,(i-1)*nr+1:i*nr) = Hess((i-1)*nr+1:i*nr,(i-1)*nr+1:i*nr)+(reg(i).*reg_temp);
    disp([all_sol(i,:) ' Regularization        = ' num2str(reg(i))]);
    fprintf(fid_log,[all_sol(i,:) ' Regularization        = %f\n'],reg(i));
  end
  clear reg_temp

  disp('Inverting Hessian');
  foo = (Hess\J'*data_diff);
  %foo = J'*minres(Hess,data_diff,1e-5,2000);

  clear J reg Hess;
  
  %******************************************************
  % Inversion complete, normalize update    
  for i = 1 : length(N)
    foo((i-1)*nr+1:i*nr) = foo((i-1)*nr+1:i*nr).*N(i);
  end
  clear nr N
  
  % Update values
  [nr,nc] = size(fwd_mesh.conc);
  foo_conc = foo(1:nr*nc);
  foo_sa = foo(nr*nc+1:nr*(nc+1));
  foo_sp = foo(nr*(nc+1)+1:end);
  
  conc_temp = conc_temp + foo_conc;
  fwd_mesh.conc = reshape(conc_temp,nr,nc);
  fwd_mesh.sa = fwd_mesh.sa + foo_sa;
  fwd_mesh.sp = fwd_mesh.sp + foo_sp;
  clear foo foo_conc foo_sa foo_sp conc_temp N
  
  %% constraining water to be less than 100%, sa and sp <3.0
  [fwd_mesh.conc, fwd_mesh.sa, fwd_mesh.sp] = ...
      constrain_val(fwd_mesh, fwd_mesh.conc, ...
		    fwd_mesh.sa, fwd_mesh.sp, ...
		    fwd_mesh.chromscattlist);
  
  %%filtering
  if (filter_n ~= 0)
    disp('Filtering');
    [fwd_mesh] = mean_filter_chromscatt(fwd_mesh, filter_n);
  end
  
  
  %**********************************************************
  % Compute absorption and scatter coefficients
  % mua From Beer's Law:
  fwd_mesh.mua = (E*fwd_mesh.conc')';
  
  % mus' From Mie Theory:
  for k =1 : nwv
    fwd_mesh.mus(:,k) = ((fwd_mesh.sa).*(wv_array(k)/1000).^(-fwd_mesh.sp));
  end
  
  % Diffusion coeff:
  fwd_mesh.kappa = (1./(3*(fwd_mesh.mus + fwd_mesh.mua)));
  
  %****************************************************************
  % Calculate Clinical Parameters
  chrom_scatt = fwd_mesh.conc;
  chrom_scatt(:,end+1) = fwd_mesh.sa;
  chrom_scatt(:,end+1) = fwd_mesh.sp;
  
  %Hbt and StO2
  %chrom_scatt(:,1) = (fwd_mesh.conc(:,1) + fwd_mesh.conc(:,2)).*1000;
  %chrom_scatt(:,2) = ((fwd_mesh.conc(:,1).*1000)./chrom_scatt(:,1))*100;
  
  %****************************************************************
  % Save solutions to file    
  
  for j = 1:n_allsol
    tmp_string = [strcat(output_fn,'_',all_sol(j,:),'.sol')];
    if (it==1)
      fid = fopen(tmp_string,'w');
    else
      fid = fopen(tmp_string,'a');
    end
    fprintf(fid,'solution %d ',it);
    fprintf(fid,'-size=%g ',length(fwd_mesh.nodes));
    fprintf(fid,'-components=1 ');
    fprintf(fid,'-type=nodal\n');
    fprintf(fid,'%f ',chrom_scatt(:,j));
    fprintf(fid,'\n');
    fclose(fid);
    clear tmp_string
  end
  clear chrom_scatt     
%  pause(2.5);
%  plotmesh(fwd_mesh);
end
% Close log file
time = toc(tmainloop);
write_log(fid_log,'Computation Time = %f\n',time,1);
write_log(fid_log,'Finished on %s ',date,1);
cl = fix(clock);
write_log(fid_log,'at %d:%d:%d\n',cl(4:6),1);
fprintf(fid_log,'Final Solution:\n');
fprintf(fid_log,'%g',it-1);
fclose(fid_log);

% END RECONSTRUCTION
%**************************************************************************




%**************************************************************************
% Selected Sub-functions

function [conc,sa,sp] = constrain_val(mesh2,conc,sa,sp,list);
% Constrain water
list = char(list);
list = list(:,1:5);
[nr,nc]=size(list);
for i = 1 : nr-2
  if strcmp(list(i,:),'Water') == 1
    index = find(conc(:,i) > 1.0);
    conc(index,i) = 1.0;
    clear index;
    index = find(conc(:,i) < 0.0);
    conc(index,i) = 0.001;
    clear index;
  else
    index = find(conc(:,i) < 0.0);
    conc(index,i) = 0.001;
    clear index;
  end
end

%%constraining scatt ampl
index = find(sa > 3.0);
sa(index) = 3.0;
clear index;   
index = find(sa < 0.0);
sa(index) = 0.1;
clear index;  

%%constraining scatt power
index = find(sp > 3.0);
sp(index) = 3.0;
clear index;     
index = find(sp < 0.0);
sp(index) = 0.1;
clear index;   
