function save_mesh(mesh,fn)

% save_mesh(mesh,fn)
%
% where mesh is the structured variable containing mesh information
% and fn is the filename to be saved
% Save mesh as seperate file of *.node, *.elem, *.param, etc
% Input mesh and filename to be saved.



% check if the user reversed the inputs
if ischar(mesh)
    temp = mesh;
    mesh = fn;
    fn = temp;
end

warning('off','MATLAB:DELETE:FileNotFound');

% saving fn.node file
mysave([fn '.node'],[mesh.bndvtx mesh.nodes]);

% saving fn.elem file
mysave([fn '.elem'],mesh.elements);

% saving fn.param file
if strcmp(mesh.type,'stnd') == 1
    mesh.kappa = 1./(3.*(mesh.mua+mesh.mus));
    data = [mesh.mua mesh.kappa mesh.ri];
elseif strcmp(mesh.type,'stnd_spn') == 1
    data = [mesh.mua mesh.mus mesh.g mesh.ri];
elseif strcmp(mesh.type,'stnd_bem') == 1
    mesh.kappa = 1./(3.*(mesh.mua+mesh.mus));
    data = [mesh.mua mesh.kappa mesh.ri];
elseif strcmp(mesh.type,'fluor') || strcmp(mesh.type,'fluor_bem')
    mesh.kappax = 1./(3.*(mesh.muax+mesh.musx));
    mesh.kappam = 1./(3.*(mesh.muam+mesh.musm));
    data = [mesh.muax mesh.kappax mesh.ri mesh.muam ...
        mesh.kappam mesh.muaf mesh.eta mesh.tau];
elseif strcmp(mesh.type,'spec') || strcmp(mesh.type,'spec_bem')
    data = [];
    for i = 1 : length(mesh.chromscattlist)
        if strcmpi(mesh.chromscattlist(i),'S-Amplitude') == 0 & ...
                strcmpi(mesh.chromscattlist(i),'S-Power') == 0
            data = [data mesh.conc(:,i)];
        end
    end
    data = [data mesh.sa mesh.sp];
end

[nrow,ncol]=size(data);
fid = fopen([fn '.param'],'w');
fprintf(fid,'%s\n',mesh.type);
if strcmp(mesh.type,'spec') || strcmp(mesh.type,'spec_bem')
    for i = 1 : length(mesh.chromscattlist)
        fprintf(fid,'%s\n',char(mesh.chromscattlist(i)));
    end
end

fmtstr = repmat('%g ',1,ncol);
fmtstr(end) = []; fmtstr = strcat(fmtstr,'\n');
fprintf(fid,fmtstr,data');
fclose(fid);
clear data

% save extinction file for spec mesh type
if strcmp(mesh.type,'spec') || strcmp(mesh.type,'spec_bem')
    
    % check if mesh.wv is the right format
    if size(mesh.wv,1)<size(mesh.wv,2)
        mesh.wv = mesh.wv';
    end
    
    data = [mesh.wv mesh.excoef];
    [nrow,ncol]=size(data);
    fid = fopen([fn '.excoef'],'w');
    for i = 1 : length(mesh.chromscattlist)
        fprintf(fid,'%s\n',char(mesh.chromscattlist(i)));
    end
    fmtstr = repmat('%g ',1,ncol);
    fmtstr(end) = []; fmtstr = strcat(fmtstr,'\n');
    fprintf(fid,fmtstr,data');
    fclose(fid);
end

% saving fn.region file
if isfield(mesh,'region') == 1
    mysave([fn '.region'],mesh.region);
else
    delete([fn '.region']);
end

% saving fn.source file
if isfield(mesh,'source') == 1
    [nrow,ncol]=size(mesh.source.coord);
    fid = fopen([fn '.source'],'w');
    if mesh.source.distributed == 1
        fprintf(fid,'%s\n','distributed');
    end
    if mesh.source.fixed == 1
        fprintf(fid,'%s\n','fixed');
    end
    if ncol == 2
        fprintf(fid, '%s\n','num    x   y   fwhm');
    elseif ncol == 3
        fprintf(fid, '%s\n','num    x   y   z  fwhm');
    end
    data = [mesh.source.num, mesh.source.coord, mesh.source.fwhm];
    [nrow,ncol] = size(data);
    fmtstr = repmat('%.12g ',1,ncol);
    fmtstr(end) = []; fmtstr = strcat(fmtstr,'\n');
    fprintf(fid,fmtstr,data');
    fclose(fid);
    clear data
else
    delete([fn '.source']);
end

% saving fn.meas file
if isfield(mesh,'meas') == 1
    [nrow,ncol]=size(mesh.meas.coord);
    fid = fopen([fn '.meas'],'w');
    if mesh.meas.fixed == 1
        fprintf(fid,'%s\n','fixed');
    end
    if ncol == 2
        fprintf(fid, '%s\n','num    x   y');
    elseif ncol == 3
        fprintf(fid, '%s\n','num    x   y   z');
    end
    data = [mesh.meas.num, mesh.meas.coord];
    [nrow,ncol] = size(data);
    fmtstr = repmat('%.12g ',1,ncol);
    fmtstr(end) = []; fmtstr = strcat(fmtstr,'\n');
    fprintf(fid,fmtstr,data');
    fclose(fid);
    clear data
else
    delete([fn '.meas']);
end

% saving fn.link file
if isfield(mesh,'link') == 1
    data = mesh.link;
    [nrow,ncol]=size(data);
    fid = fopen([fn '.link'],'w');
    fprintf(fid,'%s\n','source detector active');
    fmtstr = repmat('%g ',1,ncol);
    fmtstr(end) = []; fmtstr = strcat(fmtstr,'\n');
    fprintf(fid,fmtstr,data');
    fclose(fid);
    clear data
else
    delete([fn '.link']);
end

% saving fn.ident file
if isfield(mesh,'ident') == 1
    mysave([fn '.ident'],mesh.ident);
else
    delete([fn '.ident']);
end

warning('on','MATLAB:DELETE:FileNotFound');
