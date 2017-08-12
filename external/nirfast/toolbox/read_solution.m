function mesh=read_solution(mesh,sol_fn,it)

% mesh=read_solution(mesh,sol_fn,it)
%
% Reads solutions, for a given mesh based on the root solution file
% mesh can be either a workspace variable or a mesh filename
% sol_fn must be of the format *.sol and the correct files (as
% written by the reconstruct programs) must be present. If last
% arguement is 'it' is specified, reads solution for that given
% iteration, otherwise reads and plots final iteration



% If not a workspace variable, load mesh
if ischar(mesh)== 1
  mesh = load_mesh(mesh);
end

if strcmp(mesh.type,'stnd') || strcmp(mesh.type,'stnd_spn') || strcmp(mesh.type,'stnd_bem') % stnd
  if nargin == 2
    it = 0;
    fid = fopen([sol_fn '_mua.sol'],'r');
    while fgetl(fid) ~= -1
      it = it + 1;
    end
    it = it / 2;
    fclose(fid);
  end
  fid = fopen([sol_fn '_mua.sol'],'r');
  for i = 1 : 2*it-1
    fgetl(fid);
  end
  mesh.mua = fscanf(fid,'%g',inf);
  fclose(fid);
  fid = fopen([sol_fn '_mus.sol'],'r');
  for i = 1 : 2*it-1
    fgetl(fid);
  end
  mesh.mus = fscanf(fid,'%g',inf);
  fclose(fid);
  mesh.kappa = 1./(3.*(mesh.mua+mesh.mus));
  plotmesh(mesh);
elseif strcmp(mesh.type,'fluor') || strcmp(mesh.type,'fluor_bem')  % fluor
  if nargin == 2
    it = 0;
    fid = fopen([sol_fn '_etamuaf.sol'],'r');
    while fgetl(fid) ~= -1
      it = it + 1;
    end
    it = it / 2;
    fclose(fid);
  end
  fid = fopen([sol_fn '_etamuaf.sol'],'r');
  for i = 1 : 2*it-1
    fgetl(fid);
  end
  mesh.etamuaf = fscanf(fid,'%g',inf);
  fclose(fid);
  plotmesh(mesh);
  
elseif strcmp(mesh.type,'spec') || strcmp(mesh.type,'spec_bem')   % spec
  [nc,junk]=size(mesh.chromscattlist);
  if nargin == 2
    t = char(mesh.chromscattlist(1,1));
    it = 0;
    fid = fopen([sol_fn '_' t '.sol'],'r');
    while fgetl(fid) ~= -1
      it = it + 1;
    end
    it = it / 2;
    fclose(fid);
  end
  for j = 1 : nc-2
    t = char(mesh.chromscattlist(j,1));
    fid = fopen([sol_fn '_' t '.sol'],'r');
    for i = 1 : 2*it-1
      fgetl(fid);
    end
    mesh.conc(:,j) = fscanf(fid,'%g',inf);
    fclose(fid);
  end
  t = char(mesh.chromscattlist(j+1,1));
  fid = fopen([sol_fn '_' t '.sol'],'r');
  for i = 1 : 2*it-1
    fgetl(fid);
  end
  mesh.sa = fscanf(fid,'%g',inf);
  fclose(fid);
  t = char(mesh.chromscattlist(j+2,1));
  fid = fopen([sol_fn '_' t '.sol'],'r');
  for i = 1 : 2*it-1
    fgetl(fid);
  end
  mesh.sp = fscanf(fid,'%g',inf);
  fclose(fid);
  plotmesh(mesh)
end
