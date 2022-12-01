function export2DAE(axeshandle,filename)

[p,f]=fileparts(filename);
filename=[p f '.dae'];


fid=fopen(filename,'w');
fprintf(fid,'<?xml version="1.0"?>\n');
fprintf(fid,'<COLLADA xmlns="http://www.collada.org/2005/11/COLLADASchema" version="1.4.1">\n');
fprintf(fid,'    <asset>\n');

fprintf(fid,'        <created>2015-01-01T00:00:00</created>\n');
fprintf(fid,'        <modified>2015-01-01T00:00:00</modified>\n');
fprintf(fid,'        <up_axis>Y_UP</up_axis>\n');
fprintf(fid,'    </asset>\n');

fprintf(fid,'    <library_geometries>\n');
obj3D=findobj('parent',axeshandle,'type','patch');

for i=1:length(obj3D)
    vert=get(obj3D(i),'Vertices')';
    vnorm=-get(obj3D(i),'VertexNormals')';
    faces=get(obj3D(i),'Faces')';
    cdata=get(obj3D,'FaceVertexCData');
    cmap=colormap;
    lst=linspace(min(get(axeshandle,'clim')),max(get(axeshandle,'clim')),size(cmap,1));
    for j=1:3
        rbg(:,j)=interp1(lst,cmap(:,j),cdata);
    end
    cmap=rbg';


    fprintf(fid,'        <geometry id="shape%d-geom">\n',i-1);
    fprintf(fid,'            <mesh>\n');


    fprintf(fid,'                <source id="shape%d-geom-positions">\n',i-1);
    fprintf(fid,'                    <float_array id="shape%d-geom-positions-array" count="%d">',i-1,length(vert(:)));
    fprintf(fid,[repmat('%8.5f ', 1, length(vert(:))) '</float_array>\n'],vert(:));
    fprintf(fid,'                    <technique_common>\n');
    fprintf(fid,'                        <accessor source="#shape%d-geom-positions-array" count="%d" stride="3">\n',i-1,size(vert,2));
    fprintf(fid,'                            <param name="X" type="float"/>\n');
    fprintf(fid,'                            <param name="Y" type="float"/>\n');
    fprintf(fid,'                            <param name="Z" type="float"/>\n');
    fprintf(fid,'                        </accessor>\n');
    fprintf(fid,'                    </technique_common>\n');
    fprintf(fid,'                </source>\n');


    fprintf(fid,'                <source id="shape%d-geom-vert-colors">\n',i-1);
    fprintf(fid,'                    <float_array id="shape%d-geom-vert-colors-array" count="%d">',i-1,length(cmap(:)));
    fprintf(fid,[repmat('%8.5f ', 1, length(cmap(:))) '</float_array>\n'],cmap(:));
    fprintf(fid,'                    <technique_common>\n');
    fprintf(fid,'                        <accessor source="#shape%d-geom-vert-colors-array" count="%d" stride="3">\n',i-1,size(cmap,2));
    fprintf(fid,'                            <param name="R" type="float"/>\n');
    fprintf(fid,'                            <param name="G" type="float"/>\n');
    fprintf(fid,'                            <param name="B" type="float"/>\n');
    fprintf(fid,'                        </accessor>\n');
    fprintf(fid,'                    </technique_common>\n');
    fprintf(fid,'                </source>\n');

    fprintf(fid,'                <vertices id="shape%d-geom-vertices">\n',i-1);
    fprintf(fid,'                    <input semantic="POSITION" source="#shape%d-geom-positions"/>\n',i-1);
    fprintf(fid,'                    <input semantic="COLOR" source="#shape%d-geom-vert-colors"/>\n',i-1);
    fprintf(fid,'                </vertices>\n');

    fprintf(fid,'                <source id="shape%d-geom-normals">\n',i-1);
    fprintf(fid,'                    <float_array id="shape%d-geomnormals-array" count="%d">',i-1,length(vnorm(:)));
    fprintf(fid,[repmat('%8.5f ', 1, length(vnorm(:))) '</float_array>\n'],vnorm(:));
    fprintf(fid,'                    <technique_common>\n');
    fprintf(fid,'                        <accessor source="#shape%d-geom-normals-array" count="%d" stride="3">\n',i-1,size(vnorm,2));
    fprintf(fid,'                            <param name="X" type="float"/>\n');
    fprintf(fid,'                            <param name="Y" type="float"/>\n');
    fprintf(fid,'                            <param name="Z" type="float"/>\n');
    fprintf(fid,'                        </accessor>\n');
    fprintf(fid,'                    </technique_common>\n');
    fprintf(fid,'                </source>\n');


    fprintf(fid,'                <triangles count="%d">\n',size(faces,2));
    fprintf(fid,'                    <input semantic="VERTEX" source="#shape%d-geom-vertices" offset="0"/>\n',i-1);
    fprintf(fid,'                    <input semantic="NORMAL" source="#shape%d-geom-normals" offset="0"/>\n',i-1);
    fprintf(fid,['                    <p>' repmat('%d ',1,3*size(faces,2)) '</p>\n'],faces(:)-1);
    fprintf(fid,'                </triangles>\n');
    fprintf(fid,'            </mesh>\n');
    fprintf(fid,'        </geometry>\n');
end
fprintf(fid,'    </library_geometries>\n');

fprintf(fid,'    <library_visual_scenes>\n');
for i=1:length(obj3D)
    fprintf(fid,'        <visual_scene id="scene%d">\n',i-1);
    fprintf(fid,'            <node id="scene%d-node%d">\n',i-1,i-1);
    fprintf(fid,'                <instance_geometry url="#shape%d-geom"/>\n',i-1);
    fprintf(fid,'            </node>\n');
    fprintf(fid,'        </visual_scene>\n');
end
fprintf(fid,'    </library_visual_scenes>\n');


fprintf(fid,'    <scene>\n');
fprintf(fid,'        <instance_visual_scene url="#scene0"/>\n');
fprintf(fid,'    </scene>\n');
fprintf(fid,'</COLLADA>');


fclose(fid);







