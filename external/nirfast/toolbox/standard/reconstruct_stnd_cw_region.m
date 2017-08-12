function [fwd_mesh,pj_error] = reconstruct_stnd_cw_region(fwd_mesh,...
						  data_fn,...
						  iteration,...
						  lambda,...
						  output_fn,...
						  filter_n,...
						  region)

% [fwd_mesh,pj_error] = reconstruct_stnd_cw_region(fwd_mesh,...
%						  data_fn,...
%						  iteration,...
%						  lambda,...
%						  output_fn,...
%						  filter_n,...
%						  region)
%
% CW Reconstruction program for standard meshes using region based priori
% Needs meshes that have region labels
% Assumes that Intensity only data, and only reconstructs mua
%
% fwd_mesh is the input mesh (variable or filename)
% data_fn is the boundary data (variable or filename)
% iteration is the max number of iterations
% lambda is the initial regularization value
%    e.g. [20 20 15; 20 20 15] for a 3D mesh with 2 regions
% output_fn is the root output filename
% filter_n is the number of mean filters
% region is an array of the regions (e.g. [0 1 2])



% set modulation frequency to zero.
frequency = 0;

tic;

%****************************************
% If not a workspace variable, load mesh
if ischar(fwd_mesh)== 1
  fwd_mesh = load_mesh(fwd_mesh);
end
if ~strcmp(fwd_mesh.type,'stnd')
    errordlg('Mesh type is incorrect','NIRFAST Error');
    error('Mesh type is incorrect');
end

%*******************************************************
% read data - This is the calibrated experimental data or simulated data
anom = load_data(data_fn);
if ~isfield(anom,'paa')
    errordlg('Data not found or not properly formatted','NIRFAST Error');
    error('Data not found or not properly formatted');
end

% remove zeroed data
anom.paa(anom.link(:,3)==0,:) = [];
data_link = anom.link;

anom = anom.paa;
anom = log(anom(:,1)); %take log of amplitude
fwd_mesh.link = data_link;
%********************************************************
% Initiate projection error
pj_error = [];

%*******************************************************
% Initiate log file
fid_log = fopen([output_fn '.log'],'w');
fprintf(fid_log,'Forward Mesh   = %s\n',fwd_mesh.name);
if ischar(data_fn) ~= 0
    fprintf(fid_log,'Data File      = %s\n',data_fn);
end
fprintf(fid_log,'Initial Reg    = %d\n',lambda);
fprintf(fid_log,'Filter         = %d\n',filter_n);
fprintf(fid_log,'Output Files   = %s_mua.sol\n',output_fn);
fprintf(fid_log,'               = %s_mus.sol\n',output_fn);
fprintf(fid_log,'Initial Guess mua = %d\n',fwd_mesh.mua(1));
%*******************************************************
% This calculates the mapping matrix that reduces Jacobian from nodal
% values to regional values
disp('calculating regions');
if ~exist('region','var') || (exist('region','var') && isempty(region))
    region = unique(fwd_mesh.region);
end
K = region_mapper(fwd_mesh,region);
%*******************************************************
for it = 1 : iteration
  
  % Calculate jacobian
  [J,data]=jacobian_stnd(fwd_mesh,0);
  data.amplitude(data_link(:,3)==0,:) = [];
  J = J.complete;
  
  % reduce J into regions!
  J = J*K;

  % Read reference data
  clear ref;
  ref(:,1) = log(data.amplitude);

  data_diff = (anom-ref);

  pj_error = [pj_error sum(abs(data_diff.^2))]; 
  
  disp('---------------------------------');
  disp(['Iteration Number          = ' num2str(it)]);
  disp(['Projection error          = ' num2str(pj_error(end))]);

  fprintf(fid_log,'---------------------------------\n');
  fprintf(fid_log,'Iteration Number          = %d\n',it);
  fprintf(fid_log,'Projection error          = %f\n',pj_error(end));

  if it ~= 1
    p = (pj_error(end-1)-pj_error(end))*100/pj_error(end-1);
    disp(['Projection error change   = ' num2str(p) '%']);
    fprintf(fid_log,'Projection error change   = %f %%\n',p);
    if (p) <= 2
      disp('---------------------------------');
      disp('STOPPING CRITERIA REACHED');
      fprintf(fid_log,'---------------------------------\n');
      fprintf(fid_log,'STOPPING CRITERIA REACHED\n');
     break
    end
  end
   
   % Normalize Jacobian wrt optical values
  % This is done different to standard reconstruction, since we have region
  % based values, so we can take mean of each region.
  N = fwd_mesh.mua'*K./sum(K);
  nn = length(N);
  % Normalise by looping through each region, rather than creating a
  % diagonal matrix and then multiplying - more efficient for large meshes
  for i = 1 : nn
      J(:,i) = J(:,i).*N(i);
  end
  
  % build hessian
  [nrow,ncol]=size(J);
  Hess = zeros(ncol);
  Hess = (J'*J);
    
  % Add regularization
  % drastically reduced regularisation!
  if it ~= 1
    lambda = lambda./10;
  end

  % Regularisation is now tricky, as we want to keep it quantitative. So
  % since we are using J'*J, we will regularise for mua
  reg_mua = max(diag(Hess)).*lambda;
  
  reg = ones(ncol,1).*reg_mua;
  disp(['mua Regularization        = ' num2str(reg_mua)]);
  fprintf(fid_log,'mua Regularization        = %f\n',reg_mua);

  % Add regularisation to diagonal - looped rather than creating a matrix
  % as it is computational more efficient for large meshes
  for i = 1 : ncol
      Hess(i,i) = Hess(i,i) + reg(i);
  end
  clear reg*

  % Calculate update
  foo = Hess\(J'*data_diff);
  foo = foo.*N';
  clear nn N

  % use region mapper to unregionize!
  foo = K*foo;
  
  % Update values
  fwd_mesh.mua = fwd_mesh.mua + foo;
  fwd_mesh.kappa = 1./(3.*(fwd_mesh.mus + fwd_mesh.mua));
  
  clear foo Hess Hess_norm tmp data_diff G

  % We dont like -ve mua or mus! so if this happens, terminate
  if (any(fwd_mesh.mua<0) | any(fwd_mesh.mus<0))
    disp('---------------------------------');
    disp('-ve mua or mus calculated...not saving solution');
    fprintf(fid_log,'---------------------------------\n');
    fprintf(fid_log,'STOPPING CRITERIA REACHED\n');
    break
  end
  
  % Filtering if needed!
  if filter_n > 1
    fwd_mesh = mean_filter(fwd_mesh,abs(filter_n));
  elseif filter_n < 0
    fwd_mesh = median_filter(fwd_mesh,abs(filter_n));
  end

  if it == 1
    fid = fopen([output_fn '_mua.sol'],'w');
  else
    fid = fopen([output_fn '_mua.sol'],'a');
  end
  fprintf(fid,'solution %g ',it);
  fprintf(fid,'-size=%g ',length(fwd_mesh.nodes));
  fprintf(fid,'-components=1 ');
  fprintf(fid,'-type=nodal\n');
  fprintf(fid,'%f ',fwd_mesh.mua);
  fprintf(fid,'\n');
  fclose(fid);
  
  if it == 1
    fid = fopen([output_fn '_mus.sol'],'w');
  else
    fid = fopen([output_fn '_mus.sol'],'a');
  end
  fprintf(fid,'solution %g ',it);
  fprintf(fid,'-size=%g ',length(fwd_mesh.nodes));
  fprintf(fid,'-components=1 ');
  fprintf(fid,'-type=nodal\n');
  fprintf(fid,'%f ',fwd_mesh.mus);
  fprintf(fid,'\n');
  fclose(fid);
end

% close log file!
time = toc;
fprintf(fid_log,'Computation TimeRegularization = %f\n',time);
fclose(fid_log);




function KKK=region_mapper(mesh,region)

% created a matrix which is in effect a mapper to move 
% between nodal basis or region basis.
% h dehghani May 8 2002

%nregions = max(mesh.region)+1;
nregion = length(region);
disp(['Number of regions = ' num2str(nregion)]);
nnodes = length(mesh.nodes);

% create empty mapping matrix
K = sparse(nnodes,nregion);

% Assign mapping functions, for each node belonging to each region
for j = 1 : nregion
  K(find(mesh.region==region(j)),j) = 1;
end

% find the total number of assigned nodes
N = full(sum(sum(K)));

% Here if some node is not in region, must account for it
if N ~= length(mesh.nodes)
  KK = sparse(nnodes,nnodes-N);
  for k = 1 : length(region)
    if k == 1
      a = find(mesh.region~=region(k));
    else
      a = intersect(find(mesh.region~=region(k)),a);
    end
  end
  for i = 1 : length(a)
    KK(a(i),i) = 1;
  end
  KKK = [K KK];
else
  KKK = K;
end
