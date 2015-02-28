function plotmesh(mesh, plotflag)

% plotmesh(mesh, plotflag)
%
% Allows fast and easy viewing of mesh
% 
% mesh is the input mesh (variable or filename)
% plotflag is optional, if it is 1, the source/detectors will show




%% load mesh
if ischar(mesh)== 1
    mesh = load_mesh(mesh);
end

%% plot sources/detectors
if nargin == 1
   plotflag = 0;
end
if plotflag == 1
    figure;
    setfont(14);
    
    hold on;
    ind = find(mesh.bndvtx==1);
    if mesh.dimension == 2
        if isfield(mesh,'source') == 1
            tmp = sort(mesh.source.num);
            if length(tmp)>1
                s1 = mesh.source.coord(mesh.source.num == tmp(1),:);
                s2 = mesh.source.coord(mesh.source.num == tmp(2),:);
                plot(s1(:,1),s1(:,2),'go',s2(:,1),s2(:,2),'yo',...
                    mesh.source.coord(3:end,1),...
                    mesh.source.coord(3:end,2),'ro','LineWidth',2,'MarkerSize',8);
            else
                s1 = mesh.source.coord;
                plot(s1(:,1),s1(:,2),'go','LineWidth',2,'MarkerSize',8);
            end 
        end
        if isfield(mesh,'meas') == 1
            plot(mesh.meas.coord(:,1),...
                mesh.meas.coord(:,2),'bx','LineWidth',2,'MarkerSize',8);
        end
        plot(mesh.nodes(ind,1),mesh.nodes(ind,2),'c.');
        axis equal;
    elseif mesh.dimension == 3
        if isfield(mesh,'source') == 1
            tmp = sort(mesh.source.num);
            if length(tmp)>1
                s1 = mesh.source.coord(mesh.source.num == tmp(1),:);
                s2 = mesh.source.coord(mesh.source.num == tmp(2),:);
                plot3(s1(:,1),s1(:,2),s1(:,3),'go',s2(:,1),s2(:,2),s2(:,3),'yo',...
                    mesh.source.coord(3:end,1),...
                    mesh.source.coord(3:end,2),...
                    mesh.source.coord(3:end,3),'ro','LineWidth',2,'MarkerSize',8);
            else
                s1 = mesh.source.coord;
                plot3(s1(:,1),s1(:,2),s1(:,3),'go','LineWidth',2,'MarkerSize',8);
            end
            
        end
        if isfield(mesh,'meas') == 1
            plot3(mesh.meas.coord(:,1),...
                mesh.meas.coord(:,2),...
                mesh.meas.coord(:,3),'bx',...
                'LineWidth',2,'MarkerSize',8);
        end
        plot3(mesh.nodes(ind,1),...
            mesh.nodes(ind,2),...
            mesh.nodes(ind,3),'c.');
        axis equal;
    end
    
    if isfield(mesh,'source') && size(mesh.source.coord,1)>2
        legend('Source 1','Source 2','Sources +','Detector');
    elseif isfield(mesh,'source') && size(mesh.source.coord,1)==2
        legend('Source 1','Source 2','Detector');
    elseif isfield(mesh,'source') && size(mesh.source.coord,1)==1
        legend('Source 1','Detector');
    else
        legend('Detector');
    end
    
    xlabel('x position (mm)');
    ylabel('y position (mm)');
    zlabel('z position (mm)');
    
end



%% plot optical properties

os=computer;
if mesh.dimension == 3 ... 
        && ~strcmp(mesh.type,'stnd_bem') ... % Not BEM mesh isempty(strfind(os,'GLNX')) ... % Not Linux && ( strcmpi(os,'maci') || ~isempty(strfind(os,'PCWIN')) ) ... % Not MAC64
        && ~strcmp(mesh.type,'fluor_bem') ...
        && ~strcmp(mesh.type,'spec_bem') && ~strcmpi(os,'maci') ...
        && ~strcmp(os,'glnx86')
    
    fpath = [tempdir 'temp_nirfast.vtk'];
    disp(['Saving vtk file to: ' fpath]);
    cpwd = pwd;
    cd(fullfile(get_nirfast_root(),'toolbox/visualization/bin'))
    
    nirfast2vtk(mesh,fpath);
    if strcmpi(os,'PCWIN64')
        systemcall = ['"' which('nirviz64.exe') '" ' fpath];
    elseif strcmpi(os,'PCWIN')
        systemcall = ['"' which('nirviz.exe') '" ' fpath];
    elseif strcmpi(os,'maci64')
        nirvizcmd = which('nirviz.');
        if isempty(nirvizcmd) || ...
                isempty(regexp(nirvizcmd,'nirviz-mac', 'once'))
            nirvizcmd = GetSystemCommand('nirviz');
            if ~isempty(nirvizcmd)
                cprintf([1 0.5 0.1],...
                    ' Extracting nirviz executable!\n Please wait...\n');
                warning('off','MATLAB:dispatcher:pathWarning')
                addpath(fullfile(fileparts(nirvizcmd),'nirviz-mac'))
                warning('on','MATLAB:dispatcher:pathWarning')
                savepath
            else
                error(' Couldn''t find nirviz executable!')
            end
        end
        eval(['! chmod u+x "' nirvizcmd '"']);
        systemcall = ['DYLD_FRAMEWORK_PATH=; "' nirvizcmd '" ' fpath ...
            ' > /dev/null 2>&1 &'];
    elseif strcmpi(os,'glnxa64')
        nirvizcmd = which('nirviz.');
        if isempty(nirvizcmd) || ...
                isempty(regexp(nirvizcmd,'nirviz-linux', 'once'))
            nirvizcmd = GetSystemCommand('nirviz');
            if ~isempty(nirvizcmd)
                cprintf([1 0.5 0.1],...
                    ' Extracting nirviz executable!\n Please wait...\n');
                warning('off','MATLAB:dispatcher:pathWarning')
                addpath(fullfile(fileparts(nirvizcmd),'nirviz-linux64'))
                warning('off','MATLAB:dispatcher:pathWarning')
                savepath
            else
                error(' Couldn''t find nirviz executable!')
            end
        end
        eval(['! chmod u+x "' nirvizcmd '"']);
        systemcall = ['"' nirvizcmd '" ' fpath ' > /dev/null 2>&1 &'];
    else
        error(['OS is not supported for 3D visualization: ' computer]);
    end
    system(systemcall);
    cd(cpwd);
else

    figure;
    set(gca,'FontSize',28)

    % STANDARD
    if strcmp(mesh.type,'stnd') == 1
      subplot(1,2,1);
      plotim(mesh,mesh.mua);
      title('\mu_a','FontSize',20);
      colorbar('horiz');

      subplot(1,2,2);
      plotim(mesh,mesh.mus);
      title('\mu_s''','FontSize',20);
      colorbar('horiz');

    % STANDARD SPN
    elseif strcmp(mesh.type,'stnd_spn') == 1
      subplot(2,2,1);
      plotim(mesh,mesh.mua);
      title('\mu_a','FontSize',20);
      colorbar('horiz');

      subplot(2,2,2);
      plotim(mesh,mesh.mus);
      title('\mu_s''','FontSize',20);
      colorbar('horiz');

      subplot(2,2,3);
      plotim(mesh,mesh.g);
      title('g','FontSize',20);
      colorbar('horiz');

    % STANDARD BEM
    elseif strcmp(mesh.type,'stnd_bem') == 1
      subplot(1,2,1);
      val = zeros(size(mesh.nodes,1),1);
      for i=1:size(mesh.mua,1)
          if i==1
              ind = find(mesh.region(:,2)==0);
              ind = unique(mesh.elements(ind,:));
            val(ind) = mesh.mua(i);
          else
              ind = find(mesh.region(:,2)==i);
              ind = unique(mesh.elements(ind,:));
            val(ind) = mesh.mua(i);
          end
      end
      plotim(mesh,val);
      title('\mu_a','FontSize',20);
      colorbar('horiz');

      subplot(1,2,2);
      val = zeros(size(mesh.nodes,1),1);
      for i=1:size(mesh.mus,1)
          if i==1
              ind = find(mesh.region(:,2)==0);
              ind = unique(mesh.elements(ind,:));
            val(ind) = mesh.mus(i);
          else
              ind = find(mesh.region(:,2)==i);
              ind = unique(mesh.elements(ind,:));
            val(ind) = mesh.mus(i);
          end
      end
      plotim(mesh,val);
      title('\mu_s''','FontSize',20);
      colorbar('horiz');

    % FLUORESCENCE BEM
    elseif strcmp(mesh.type,'fluor_bem') == 1 
      subplot(3,2,1);
      val = zeros(size(mesh.nodes,1),1);
      for i=1:size(mesh.muax,1)
          if i==1
              ind = find(mesh.region(:,2)==0);
              ind = unique(mesh.elements(ind,:));
            val(ind) = mesh.muax(i);
          else
              ind = find(mesh.region(:,2)==i);
              ind = unique(mesh.elements(ind,:));
            val(ind) = mesh.muax(i);
          end
      end
      plotim(mesh,val);
      title('\mu_{ax}','FontSize',10);
      colorbar('horiz');

      subplot(3,2,2);
      val = zeros(size(mesh.nodes,1),1);
      for i=1:size(mesh.musx,1)
          if i==1
              ind = find(mesh.region(:,2)==0);
              ind = unique(mesh.elements(ind,:));
            val(ind) = mesh.musx(i);
          else
              ind = find(mesh.region(:,2)==i);
              ind = unique(mesh.elements(ind,:));
            val(ind) = mesh.musx(i);
          end
      end
      plotim(mesh,val);
      title('\mu_{sx}','FontSize',10);
      colorbar('horiz');

      subplot(3,2,3);
      val = zeros(size(mesh.nodes,1),1);
      for i=1:size(mesh.muam,1)
          if i==1
              ind = find(mesh.region(:,2)==0);
              ind = unique(mesh.elements(ind,:));
            val(ind) = mesh.muam(i);
          else
              ind = find(mesh.region(:,2)==i);
              ind = unique(mesh.elements(ind,:));
            val(ind) = mesh.muam(i);
          end
      end
      plotim(mesh,val);
      title('\mu_{am}','FontSize',10);
      colorbar('horiz');

      subplot(3,2,4);
      val = zeros(size(mesh.nodes,1),1);
      for i=1:size(mesh.musm,1)
          if i==1
              ind = find(mesh.region(:,2)==0);
              ind = unique(mesh.elements(ind,:));
            val(ind) = mesh.musm(i);
          else
              ind = find(mesh.region(:,2)==i);
              ind = unique(mesh.elements(ind,:));
            val(ind) = mesh.musm(i);
          end
      end
      plotim(mesh,val);
      title('\mu_{sm}','FontSize',10);
      colorbar('horiz');

      if ~isfield(mesh,'etamuaf')
          mesh.etamuaf = mesh.muaf.*mesh.eta;
      end

      subplot(3,2,5);
      val = zeros(size(mesh.nodes,1),1);
      for i=1:size(mesh.etamuaf,1)
          if i==1
              ind = find(mesh.region(:,2)==0);
              ind = unique(mesh.elements(ind,:));
            val(ind) = mesh.etamuaf(i);
          else
              ind = find(mesh.region(:,2)==i);
              ind = unique(mesh.elements(ind,:));
            val(ind) = mesh.etamuaf(i);
          end
      end
      plotim(mesh,val);
      title('\eta\mu_{fl}','FontSize',10);
      colorbar('horiz');

      subplot(3,2,6);
      val = zeros(size(mesh.nodes,1),1);
      for i=1:size(mesh.tau,1)
          if i==1
              ind = find(mesh.region(:,2)==0);
              ind = unique(mesh.elements(ind,:));
            val(ind) = mesh.tau(i);
          else
              ind = find(mesh.region(:,2)==i);
              ind = unique(mesh.elements(ind,:));
            val(ind) = mesh.tau(i);
          end
      end
      plotim(mesh,val);
      title('\tau','FontSize',10);
      colorbar('horiz');

    % FLUORESCENCE
    elseif strcmp(mesh.type,'fluor') == 1
      subplot(3,2,1);
      plotim(mesh,mesh.muax);
      title('\mu_{ax}','FontSize',10);
      colorbar;

      subplot(3,2,2);
      plotim(mesh,mesh.musx);
      title('\mu_{sx}''','FontSize',10);
      colorbar;

      subplot(3,2,3);
      plotim(mesh,mesh.muam);
      title('\mu_{am}','FontSize',10);
      colorbar;

      subplot(3,2,4);
      plotim(mesh,mesh.musm);
      title('\mu_{sm}''','FontSize',10);
      colorbar;

      subplot(3,2,5);
      if isfield(mesh,'etamuaf') == 1
          plotim(mesh,mesh.etamuaf);
      else
          plotim(mesh,mesh.muaf.*mesh.eta);
      end
      title('\eta\mu_{fl}','FontSize',10);
      colorbar;

      subplot(3,2,6);
      plotim(mesh,mesh.tau);
      title('\tau','FontSize',10);
      colorbar;

    % SPECTRAL BEM
    elseif strcmp(mesh.type,'spec_bem') == 1
      [nc,junk]=size(mesh.chromscattlist);

      if isfield(mesh,'etamuaf')
        n = ceil((nc-2)/2)+2;
      else
        n = ceil((nc-2)/2)+1;
      end
      k = 0;
      for i = 1 : nc-2
        k = k + 1;
        subplot(n,2,k);

          tmpv = mesh.conc(:,i);
          val = zeros(size(mesh.nodes,1),1);
          for iii=1:size(tmpv,1)
              if iii==1
                  ind = find(mesh.region(:,2)==0);
                  ind = unique(mesh.elements(ind,:));
                val(ind) = tmpv(iii);
              else
                  ind = find(mesh.region(:,2)==iii);
                  ind = unique(mesh.elements(ind,:));
                val(ind) = tmpv(iii);
              end
          end

        plotim(mesh,val);
        t = char(mesh.chromscattlist(i,1));
        title(t,'FontSize',10);
        colorbar;
      end
      subplot(n,2,k+1);
      val = zeros(size(mesh.nodes,1),1);
          for i=1:size(mesh.sa,1)
              if i==1
                  ind = find(mesh.region(:,2)==0);
                  ind = unique(mesh.elements(ind,:));
                val(ind) = mesh.sa(i);
              else
                  ind = find(mesh.region(:,2)==i);
                  ind = unique(mesh.elements(ind,:));
                val(ind) = mesh.sa(i);
              end
          end
        plotim(mesh,val);
      title('Scatter Amplitude','FontSize',10);
      colorbar;
      subplot(n,2,k+2);
      val = zeros(size(mesh.nodes,1),1);
          for i=1:size(mesh.sp,1)
              if i==1
                  ind = find(mesh.region(:,2)==0);
                  ind = unique(mesh.elements(ind,:));
                val(ind) = mesh.sp(i);
              else
                  ind = find(mesh.region(:,2)==i);
                  ind = unique(mesh.elements(ind,:));
                val(ind) = mesh.sp(i);
              end
          end
        plotim(mesh,val);
      title('Scatter Power','FontSize',10);
      colorbar;
      if isfield(mesh,'etamuaf')
        subplot(n,2,k+3);
        val = zeros(size(mesh.nodes,1),1);
          for i=1:size(mesh.etamuaf,1)
              if i==1
                  ind = find(mesh.region(:,2)==0);
                  ind = unique(mesh.elements(ind,:));
                val(ind) = mesh.etamuaf(i);
              else
                  ind = find(mesh.region(:,2)==i);
                  ind = unique(mesh.elements(ind,:));
                val(ind) = mesh.etamuaf(i);
              end
          end
        plotim(mesh,val);
        title('etamuaf','FontSize',10);
        colorbar;
      end

    % SPECTRAL
    elseif strcmp(mesh.type,'spec') == 1
      [nc,junk]=size(mesh.chromscattlist);

      if isfield(mesh,'etamuaf')
        n = ceil((nc-2)/2)+2;
      else
        n = ceil((nc-2)/2)+1;
      end
      k = 0;
      for i = 1 : nc-2
        k = k + 1;
        subplot(n,2,k);
        plotim(mesh,mesh.conc(:,i));
        t = char(mesh.chromscattlist(i,1));
        title(t,'FontSize',10);
        colorbar;
      end
      subplot(n,2,k+1);
      plotim(mesh,mesh.sa);
      title('Scatter Amplitude','FontSize',10);
      colorbar;
      subplot(n,2,k+2);
      plotim(mesh,mesh.sp);
      title('Scatter Power','FontSize',10);
      colorbar;
      if isfield(mesh,'etamuaf')
        subplot(n,2,k+3);
        plotim(mesh,mesh.etamuaf);
        title('etamuaf','FontSize',10);
        colorbar;
      end
    else
        errordlg('Mesh type not supported','NIRFAST Error');
        error('Mesh type not supported');
    end

end

end


%% plot image function
function plotim(mesh,val)
if mesh.dimension == 3 && ...
        (strcmp(mesh.type,'stnd_bem') || strcmp(mesh.type,'spec_bem') || strcmp(mesh.type,'fluor_bem'))
    h = trisurf(mesh.elements,...
	    mesh.nodes(:,1),...
	    mesh.nodes(:,2),...
	    mesh.nodes(:,3),...
	    val,'FaceAlpha',0.5);
else
    h = trisurf(mesh.elements,...
	    mesh.nodes(:,1),...
	    mesh.nodes(:,2),...
	    mesh.nodes(:,3),...
	    val);
end
shading interp;
view(2);
axis equal; 
axis off;
colormap hot;

end

%% set font function
function setfont(s)

title([],'fontsize',s);
xlabel([],'fontsize',s);
xlabel([],'fontsize',s);
set(gca,'fontsize',s);
        
end
