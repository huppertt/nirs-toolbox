function probe2xml(probe,filename)

fid = fopen(filename,'w');
fprintf(fid,'<?xml version="1.0" encoding="UTF-8" ?>\n');
fprintf(fid,'<probe>\n');

wav=unique(probe.link.type);
for i=1:length(wav)
    fprintf(fid,'\t<wavelength> %d </wavelength>\n',wav(i));
end

fprintf(fid,'\n');

for i=1:size(probe.srcPos,1);
    fprintf(fid,'\t<srcpos> <index> %d </index><x> %g </x> <y> %g </y> </srcpos>\n',...
        round(i),probe.srcPos(i,1),probe.srcPos(i,2));
end

if(isa(probe,'nirs.core.Probe1020'))
    fprintf(fid,'\n');

    p2=probe.swap_reg;
    
    for i=1:size(p2.srcPos,1);
        fprintf(fid,'\t<srcpos3D> <index> %d </index> <x> %g </x> <y> %g </y> <z> %g </z> </srcpos3D>\n',...
            round(i),p2.srcPos(i,1),p2.srcPos(i,2),p2.srcPos(i,3));
    end
end
fprintf(fid,'\n');


for i=1:size(probe.detPos,1);
    fprintf(fid,'\t<detpos> <index> %d </index> <x> %g </x> <y> %g </y> </detpos>\n',...
        round(i),probe.detPos(i,1),probe.detPos(i,2));
end

if(isa(probe,'nirs.core.Probe1020'))
    fprintf(fid,'\n');

    p2=probe.swap_reg;
    
    for i=1:size(p2.detPos,1);
        fprintf(fid,'\t<detpos3D> <index> %d </index> <x> %g </x> <y> %g </y> <z> %g </z> </detpos3D>\n',...
            round(i),p2.detPos(i,1),p2.detPos(i,2),p2.detPos(i,3));
    end
end
fprintf(fid,'\n');


ml=round(unique([probe.link.source probe.link.detector],'rows'));
for i=1:size(ml,1)
      fprintf(fid,'\t<ml> <src> %d </src> <det> %d </det> </ml>\n',ml(i,1),ml(i,2));
end

fprintf(fid,'</probe>');
fclose(fid);