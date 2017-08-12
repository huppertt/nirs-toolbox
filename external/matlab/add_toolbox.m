
name=which(name);
name2=name(length('/Applications/MATLAB_R2015a.app/toolbox/')+1:end);
p2='/Users/admin/Desktop/nirs-toolbox/external/matlab';

pp=fileparts(name2);
system(['mkdir -p ' fullfile(p2,'toolbox',pp)]);
copyfile(name,fullfile(p2,'toolbox',name2));