function export2DAE(axeshandle,filename)

[p,f]=fileparts(filename);

obj3D=findobj('parent',axeshandle,'type','patch');

for i=1:length(obj3D)
    vert=get(obj3D(i),'Vertices');
    vnorm=get(obj3D(i),'VertexNormals');
    faces=get(obj3D(i),'Faces');
    cdata=get(obj3D(i),'CData');
    if(isempty(cdata));
        cdata=ones(size(faces,1),1)*get(obj3D(i),'FaceColor');
    else
        fvcdata=get(obj3D(i),'FaceVertexCData');
        cmap=colormap;
        lst=linspace(min(get(axeshandle,'clim')),max(get(axeshandle,'clim')),size(cmap,1));
        for j=1:3
            rbg(:,1,j)=interp1(lst,cmap(:,1),fvcdata(cdata(j,:))');
            rbg(:,2,j)=interp1(lst,cmap(:,2),fvcdata(cdata(j,:))');
            rbg(:,3,j)=interp1(lst,cmap(:,3),fvcdata(cdata(j,:))');
        end
        cdata=mean(rbg,3);
    end
    png=zeros(ceil(sqrt(size(faces,1)))*ceil(sqrt(size(faces,1))),3);
    png(1:length(cdata),:)=cdata;

    png=reshape(png,[ceil(sqrt(size(faces,1))), ceil(sqrt(size(faces,1))),3]);
    [T,S]=ind2sub([ceil(sqrt(size(faces,1))), ceil(sqrt(size(faces,1)))],1:length(cdata));
    ST(:,1)=S/ceil(sqrt(size(faces,1))); ST(:,2)=T/ceil(sqrt(size(faces,1)));
end


vert=vert';
vnorm=-vnorm';
faces=faces';
ST=ST';

filename=[p f '.dae'];
pngfile=[p f '.png'];

imwrite(png,pngfile);


fid=fopen(filename,'w');
fprintf(fid,'<?xml version="1.0" encoding="UTF-8"?>\n');
fprintf(fid,'<COLLADA xmlns="http://www.collada.org/2005/11/COLLADASchema" version="1.4.1">\n');
fprintf(fid,'<asset>\n');
%  <contributor>
%   <authoring_tool>SceneKit Collada Exporter v1.0</authoring_tool>
%  </contributor>
%  <created>2014-07-02T15:29:44Z</created>
%  <modified>2014-07-02T15:29:44Z</modified>
fprintf(fid,'\t<unit meter="0.010000"/>\n');
fprintf(fid,'\t<up_axis>Y_UP</up_axis>\n');
fprintf(fid,'</asset>\n');
fprintf(fid,'<library_images>\n');
fprintf(fid,'\t<image id="image1">\n');
fprintf(fid,'\t\t<init_from>%s</init_from>\n',pngfile);
fprintf(fid,'\t</image>\n');
fprintf(fid,'</library_images>\n');

fprintf(fid,'<library_materials>\n');
fprintf(fid,'\t<material id="lambert1" name="lambert1">\n');
fprintf(fid,'\t\t<instance_effect url="#effect_lambert1"/>\n');
fprintf(fid,'\t</material>\n');
fprintf(fid,'</library_materials>\n');

fprintf(fid,'<library_effects>\n');
fprintf(fid,'\t<effect id="effect_lambert1">\n');
fprintf(fid,'\t\t<profile_COMMON>\n');
fprintf(fid,'\t\t\t<newparam sid="ID2_image1_surface">\n');
fprintf(fid,'\t\t\t\t<surface type="2D">\n');
fprintf(fid,'\t\t\t\t\t<init_from>image1</init_from>\n');
fprintf(fid,'\t\t\t\t</surface>\n');
fprintf(fid,'\t\t\t</newparam>\n');
fprintf(fid,'\t\t\t<newparam sid="ID2_image1">\n');
fprintf(fid,'\t\t\t\t<sampler2D>\n');
fprintf(fid,'\t\t\t\t\t<source>ID2_image1_surface</source>\n');
fprintf(fid,'\t\t\t\t\t<wrap_s>CLAMP</wrap_s>\n');
fprintf(fid,'\t\t\t\t\t<wrap_t>CLAMP</wrap_t>\n');
fprintf(fid,'\t\t\t\t\t<minfilter>LINEAR</minfilter>\n');
fprintf(fid,'\t\t\t\t\t<magfilter>LINEAR</magfilter>\n');
fprintf(fid,'\t\t\t\t\t<mipfilter>LINEAR</mipfilter>\n');
fprintf(fid,'\t\t\t\t</sampler2D>\n');
fprintf(fid,'\t\t\t</newparam>\n');
fprintf(fid,'\t\t\t<technique sid="common">\n');
fprintf(fid,'\t\t\t\t<blinn>\n');
fprintf(fid,'\t\t\t\t\t<ambient>\n');
fprintf(fid,'\t\t\t\t\t\t<color>1 1 1 1</color>\n');
fprintf(fid,'\t\t\t\t\t</ambient>\n');
fprintf(fid,'\t\t\t\t\t<diffuse>\n');
fprintf(fid,'\t\t\t\t\t\t<texture texture="ID2_image1" texcoord="CHANNEL2"/>\n');
fprintf(fid,'\t\t\t\t\t</diffuse>\n');
fprintf(fid,'\t\t\t\t\t<specular>\n');
fprintf(fid,'\t\t\t\t\t\t<color>0.496564 0.496564 0.496564 1</color>\n');
fprintf(fid,'\t\t\t\t\t</specular>\n');
fprintf(fid,'\t\t\t\t\t<shininess>\n');
fprintf(fid,'\t\t\t\t\t\t<float>0.022516</float>\n');
fprintf(fid,'\t\t\t\t\t</shininess>\n');
fprintf(fid,'\t\t\t\t\t<reflective>\n');
fprintf(fid,'\t\t\t\t\t\t<color>0 0 0 1</color>\n');
fprintf(fid,'\t\t\t\t\t</reflective>\n');
fprintf(fid,'\t\t\t\t\t<transparent opaque="A_ONE">\n');
fprintf(fid,'\t\t\t\t\t\t<color>0.995918 1 1 1</color>\n');
fprintf(fid,'\t\t\t\t\t</transparent>\n');
fprintf(fid,'\t\t\t\t\t<transparency>\n');
fprintf(fid,'\t\t\t\t\t\t<float>1</float>\n');
fprintf(fid,'\t\t\t\t\t</transparency>\n');
fprintf(fid,'\t\t\t\t\t<index_of_refraction>\n');
fprintf(fid,'\t\t\t\t\t\t<float>1</float>\n');
fprintf(fid,'\t\t\t\t\t</index_of_refraction>\n');
fprintf(fid,'\t\t\t\t</blinn>\n');
fprintf(fid,'\t\t\t</technique>\n');
fprintf(fid,'\t\t</profile_COMMON>\n');
fprintf(fid,'\t\t<extra>\n');
fprintf(fid,'\t\t\t<technique profile="SceneKit">\n');
fprintf(fid,'\t\t\t\t<litPerPixel>1</litPerPixel>\n');
fprintf(fid,'\t\t\t\t<ambient_diffuse_lock>1</ambient_diffuse_lock>\n');
fprintf(fid,'\t\t\t\t<intensities>\n');
fprintf(fid,'\t\t\t\t\t<emission>\n');
fprintf(fid,'\t\t\t\t\t\t<float>0.5</float>\n');
fprintf(fid,'\t\t\t\t\t</emission>\n');
fprintf(fid,'\t\t\t\t</intensities>\n');
fprintf(fid,'\t\t\t</technique>\n');
fprintf(fid,'\t\t</extra>\n');
fprintf(fid,'\t</effect>\n');
fprintf(fid,'</library_effects>\n');

fprintf(fid,'<library_geometries>\n');
fprintf(fid,'\t<geometry id="Scrap_MeshShape" name="Scrap_MeshShape">\n');
fprintf(fid,'\t\t<extra>\n');
fprintf(fid,'\t\t\t<technique profile="SceneKit">\n');
fprintf(fid,'\t\t\t\t<double_sided>1</double_sided>\n');
fprintf(fid,'\t\t\t</technique>\n');
fprintf(fid,'\t\t</extra>\n');


fprintf(fid,'\t\t<mesh>\n');
fprintf(fid,'\t\t\t<source id="Scrap_MeshShape-positions">\n');
fprintf(fid,['\t\t\t\t<float_array id="ID3-array" count="%d">' repmat(' %8.5f ', 1, 3*size(vert,2)) '</float_array>\n'],3*size(vert,2),vert(:));
fprintf(fid,'\t\t\t\t<technique_common>\n');
fprintf(fid,'\t\t\t\t\t<accessor source="#ID3-array" count="%d" stride="3">\n',size(vert,2));
fprintf(fid,'\t\t\t\t\t\t<param name="X" type="float"/>\n');
fprintf(fid,'\t\t\t\t\t\t<param name="Y" type="float"/>\n');
fprintf(fid,'\t\t\t\t\t\t<param name="Z" type="float"/>\n');
fprintf(fid,'\t\t\t\t\t</accessor>\n');
fprintf(fid,'\t\t\t\t</technique_common>\n');
fprintf(fid,'\t\t\t</source>\n');

fprintf(fid,'\t\t\t<source id="Scrap_MeshShape-normals">\n');
fprintf(fid,['\t\t\t\t<float_array id="ID4-array" count="%d">' repmat(' %8.5f ', 1, 3*size(vnorm,2)) '</float_array>\n'],3*size(vnorm,2),vnorm(:));      
fprintf(fid,'\t\t\t\t<technique_common>\n');
fprintf(fid,'\t\t\t\t\t<accessor source="#ID4-array" count="%d" stride="3">\n',size(vnorm,2));
fprintf(fid,'\t\t\t\t\t\t<param name="X" type="float"/>\n');
fprintf(fid,'\t\t\t\t\t\t<param name="Y" type="float"/>\n');
fprintf(fid,'\t\t\t\t\t\t<param name="Z" type="float"/>\n');
fprintf(fid,'\t\t\t\t\t</accessor>\n');
fprintf(fid,'\t\t\t\t</technique_common>\n');
fprintf(fid,'\t\t\t</source>\n');

fprintf(fid,'\t\t\t<source id="Scrap_MeshShape-map1">\n');
fprintf(fid,['\t\t\t\t<float_array id="ID5-array" count="%d">' repmat(' %8.5f ', 1, 2*size(ST,2)) '</float_array>\n'],2*size(ST,2),ST(:));      
fprintf(fid,'\t\t\t\t<technique_common>\n');
fprintf(fid,'\t\t\t\t\t<accessor source="#ID4-array" count="%d" stride="2">\n',size(ST,2));
fprintf(fid,'\t\t\t\t\t\t<param name="S" type="float"/>\n');
fprintf(fid,'\t\t\t\t\t\t<param name="T" type="float"/>\n');
fprintf(fid,'\t\t\t\t\t</accessor>\n');
fprintf(fid,'\t\t\t\t</technique_common>\n');
fprintf(fid,'\t\t\t</source>\n');


fprintf(fid,'\t\t\t<vertices id="Scrap_MeshShape-positions-vertices">\n');
fprintf(fid,'\t\t\t\t<input semantic="POSITION" source="#Scrap_MeshShape-positions"/>\n');
fprintf(fid,'\t\t\t</vertices>\n');

fprintf(fid,'\t\t\t<triangles count="%d" material="geometryElement6">\n',size(faces,2));
fprintf(fid,'\t\t\t\t<input semantic="VERTEX" offset="0" source="#Scrap_MeshShape-positions-vertices"/>\n');
fprintf(fid,'\t\t\t\t<input semantic="NORMAL" offset="0" source="#Scrap_MeshShape-normals"/>\n');
fprintf(fid,'\t\t\t\t<input semantic="TEXCOORD" offset="0" source="#Scrap_MeshShape-map1" set="1"/>\n');
fprintf(fid,['\t\t\t\t<p>' repmat(' %d ',1,3*size(faces,2)) '</p>\n'],faces(:)-1);
fprintf(fid,'\t\t\t</triangles>\n');
fprintf(fid,'\t\t</mesh>\n');
fprintf(fid,'\t</geometry>\n');
fprintf(fid,'</library_geometries>\n');

fprintf(fid,'<library_visual_scenes>\n');
fprintf(fid,'\t<visual_scene id="scene7">\n');
fprintf(fid,'\t\t<node id="ship" name="ship">\n');
fprintf(fid,'\t\t\t<node id="shipMesh" name="shipMesh">\n');
fprintf(fid,'\t\t\t\t<matrix>0.2 0 0 0 0 0.2 0 -1.187411 0 0 0.2 0.7991223 0 0 0 1 </matrix>\n');
fprintf(fid,'\t\t\t\t<instance_geometry url="#Scrap_MeshShape">\n');
fprintf(fid,'\t\t\t\t\t<bind_material>\n');
fprintf(fid,'\t\t\t\t\t\t<technique_common>\n');
fprintf(fid,'\t\t\t\t\t\t\t<instance_material symbol="geometryElement6" target="#lambert1">\n');
fprintf(fid,'\t\t\t\t\t\t\t\t<bind_vertex_input semantic="CHANNEL2" input_semantic="TEXCOORD" input_set="1"/>\n');
fprintf(fid,'\t\t\t\t\t\t\t</instance_material>\n');
fprintf(fid,'\t\t\t\t\t\t</technique_common>\n');
fprintf(fid,'\t\t\t\t\t</bind_material>\n');
fprintf(fid,'\t\t\t\t</instance_geometry>\n');
fprintf(fid,'\t\t\t\t<node id="emitter" name="emitter">\n');
fprintf(fid,'\t\t\t\t\t<matrix>-0.9998822 0.01057715 -0.0111275 0.009105275 0.01387102 0.31177 -0.9500563 5.406124 -0.006579665 -0.9500988 -0.3118801 -22.78586 0 0 0 1 </matrix>\n');
fprintf(fid,'\t\t\t\t</node>\n');
fprintf(fid,'\t\t\t</node>\n');
fprintf(fid,'\t\t</node>\n');
fprintf(fid,'\t</visual_scene>\n');
fprintf(fid,'</library_visual_scenes>\n');
fprintf(fid,'<scene>\n');
fprintf(fid,'\t<instance_visual_scene url="#scene7"/>\n');
fprintf(fid,'</scene>\n');
fprintf(fid,'</COLLADA>');

fclose(fid);


    
    



