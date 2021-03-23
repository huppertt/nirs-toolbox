function [aux,t,hdr]=loadBrainVisionAux(filename)

[hdr.fs, hdr.label,hdr.meta]=bva_readheader(filename);
d=bva_loadeeg(filename);
aux=d; %(33:end,:);
t=[0:size(aux,2)-1]./hdr.fs;

return