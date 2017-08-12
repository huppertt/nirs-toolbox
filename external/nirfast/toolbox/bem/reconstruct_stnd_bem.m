function [fwd_mesh,pj_error] = reconstruct_stnd_bem(fwd_mesh,...
                                                frequency,...
                                                data_fn,...
                                                iteration,...
                                                lambda,...
                                                output_fn,...
                                                filter_n)

% [fwd_mesh,pj_error] = reconstruct_stnd_bem(fwd_mesh,...
%                                        frequency,...
%                                        data_fn,...
%                                        iteration,...
%                                        lambda,...
%                                        output_fn,...
%                                        filter_n)
%                                            
% Reconstruction program for standard bem meshes
%
% fwd_mesh is the input mesh (variable or filename)
% frequency is the modulation frequency (MHz)
% data_fn is the boundary data (variable or filename)
% iteration is the max number of iterations
% lambda is the initial regularization value
% output_fn is the root output filename
% filter_n is the number of mean filters



tmainloop = tic;

if frequency < 0
    errordlg('Frequency must be nonnegative','NIRFAST Error');
    error('Frequency must be nonnegative');
end

%****************************************
% If not a workspace variable, load mesh
if ischar(fwd_mesh)== 1
  fwd_mesh = load_mesh(fwd_mesh);
end

if ~strcmp(fwd_mesh.type,'stnd_bem')
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
anom(:,1) = log(anom(:,1)); %take log of amplitude
anom(:,2) = anom(:,2)/180.0*pi; % phase is in radians and not degrees
anom(anom(:,2)<0,2) = anom(anom(:,2)<0,2) + (2*pi);
anom(anom(:,2)>(2*pi),2) = anom(anom(:,2)>(2*pi),2) - (2*pi);
anom = reshape(anom',size(anom,1)*2,1); 

fwd_mesh.link = data_link;
clear data
%******************************************************* 

% Initiate projection error
pj_error=zeros(1,iteration);

% Initiate log file
fid_log = fopen([output_fn '.log'],'w');
cl = fix(clock);
write_log(fid_log,'Started on %s ',date,1);
write_log(fid_log,'at %d:%d:%d\n',cl(4:6),1);
write_log(fid_log,'Forward Mesh   = %s\n',fwd_mesh.name);
write_log(fid_log,'Frequency      = %f MHz\n',frequency);
if ischar(data_fn) ~= 0
    write_log(fid_log,'Data File      = %s\n',data_fn);
end
if isstruct(lambda)
    write_log(fid_log,'Initial Regularization  = %d\n',lambda.value);
else
    write_log(fid_log,'Initial Regularization  = %d\n',lambda);
end
write_log(fid_log,'Filter         = %d\n',filter_n);
write_log(fid_log,'Output Files   = %s_mua.sol\n',output_fn);
write_log(fid_log,'               = %s_mus.sol\n',output_fn);
fprintf(fid_log,'Initial Guess mua = %d\n',fwd_mesh.mua(1));
fprintf(fid_log,'Initial Guess mus = %d\n',fwd_mesh.mus(1));


% start non-linear itertaion image reconstruction part
for it = 1 : iteration
  
  % Calculate jacobian
  [J,data]=jacobian_stnd_bem(fwd_mesh,frequency);
  data.paa(data_link(:,3)==0,:) = [];

  % Read reference data calculated by initial -current- guess
  clear ref;
  ref = data.paa;
  ref = reshape(ref',size(ref,1)*2,1);
  
  if size(anom,1) ~= size(ref,1)
      errordlg('Data size is incorrect','NIRFAST Error');
      error('Data size is incorrect');
  end
  data_diff = (anom-ref);

  % PJ error
  pj_error(it) = sum(abs(data_diff.^2));
 
  msg={'---------------------------------\n';...
       'Iteration Number          = %d\n';...
       'Projection error          = %f\n'};
  msgdata{1}=[]; msgdata{2}=it;msgdata{3}=pj_error(it);
  write_log(fid_log,msg,msgdata,1);
  
  if it ~= 1
    p = (pj_error(it-1)-pj_error(it))*100/pj_error(it-1);
    write_log(fid_log,'Projection error change   = %f %%\n',p,1);
    if p <= 0.5 || (pj_error(it) < (10^-18)) % stopping criteria is currently set at 2% decrease change
      write_log(fid_log,{'---------------------------------\n';...
                         'STOPPING CRITERIA REACHED\n'},{[];[]},1);
     break
    end
  end
  
  % Normalize Jacobian
  J = J*diag([fwd_mesh.mua;fwd_mesh.kappa]);
  
  % Add regularization, which decreases at each iteration
  if it ~= 1
    lambda = lambda/10^0.25;
  end
  
  % build hessian
  Hess = (J'*J);

  reg = lambda*max(diag(Hess));
  write_log(fid_log,'Regularization            = %f\n',reg,1);
  Hess = Hess + reg*eye(length(Hess));
  data_diff = J'*data_diff;
 
  % Calculate update
  foo = (Hess\data_diff);
  delta_mua = foo(1:end/2).*fwd_mesh.mua;
  delta_kappa = foo(end/2+1:end).*fwd_mesh.kappa;

  % Update values
  fwd_mesh.mua = fwd_mesh.mua + delta_mua;
  fwd_mesh.kappa = fwd_mesh.kappa + delta_kappa;
  fwd_mesh.mus = (1./(3.*fwd_mesh.kappa)) - fwd_mesh.mua;
  
  clear foo Hess Hess_norm tmp data_diff G
  
  % We dont like -ve mua or mus! so if this happens, terminate
  if (any(fwd_mesh.mua<0) || any(fwd_mesh.mus<0))
    write_log(fid_log,{'Negative mua or mus calculated...not saving solution';...
                       '---------------------------------\n';...
                       'STOPPING CRITERIA REACHED\n'},{[];[];[]},1);
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
time = toc(tmainloop);
write_log(fid_log,'Computation Time = %f\n',time,1);
write_log(fid_log,'Finished on %s ',date,1);
cl = fix(clock);
write_log(fid_log,'at %d:%d:%d\n',cl(4:6),1);
fclose(fid_log);



function write_log(fid,msg,msgdata,printflag)
% writes the info in msgdata formatted as msg into file whose id is 'fid'
% if 4th argument is present it will also print on screen.
if nargin<=3
    printflag=0;
end

if ~iscell(msg)
    mymsg={msg};
    mymsgdata={msgdata};
else
    mymsg=msg;
    mymsgdata=msgdata;
end

for i=1:length(mymsg)
    fprintf(fid,mymsg{i},mymsgdata{i});
    if printflag
        fprintf(mymsg{i},mymsgdata{i});
    end
end
