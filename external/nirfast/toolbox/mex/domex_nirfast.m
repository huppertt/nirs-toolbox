function domex_nirfast()
% Routine to compile the meshing mex files in a batch.
% Should only be run on a machine that has mex set up properly.
% It should be called whenever meshing mex C/C++ files are updated.
% A normal user is unlikely to run this script!
% Written by Hamid Ghadyani June, 2010

mypwd = pwd;
% Find the folder where meshing mex files reside
meshing_mex = fileparts(which('nirfast'));
cd(meshing_mex);
cd('../toolbox/mex/meshing_mex')
mexcommand = 'mex -v ';
mexoptions = ' '; % use this for debugging
os=computer;
if isunix
    if ~isempty(strfind(os,'MAC')) % Mac
        eval([mexcommand mexoptions ' -DOSX -I./meshlib tag_checkerboard3d_mex.cpp meshlib/CStopWatch.cpp'])
        eval([mexcommand mexoptions ' -I./meshlib             surface_relations_mex.cpp isinvolume_randRay.cpp meshlib/geomath.cpp meshlib/vector.cpp meshlib/polyhedron2BSP.cpp meshlib/CPoint.cpp meshlib/CVector.cpp meshlib/Plane3D.cpp meshlib/BSPNode.cpp meshlib/MeshIO.cpp meshlib/FileOperation.cpp meshlib/CPolygon.cpp meshlib/predicates.cpp meshlib/CStopWatch.cpp'])
        eval([mexcommand mexoptions ' pnpoly_mex.cpp'])
        eval([mexcommand mexoptions '        -I./meshlib orient_surface_mex.cpp meshlib/vector.cpp meshlib/geomath.cpp isinvolume_randRay.cpp'])
        eval([mexcommand mexoptions ' -I./meshlib involume_mex.cpp isinvolume_randRay.cpp meshlib/geomath.cpp meshlib/vector.cpp'])
        eval([mexcommand mexoptions '         -I./meshlib intersect_ray_shell_mex.cpp ./meshlib/vector.cpp ./meshlib/geomath.cpp'])
        eval([mexcommand mexoptions ' extract_ind_regions_mex.cpp'])
        eval([mexcommand mexoptions ' -I./meshlib         PointInPolyhedron_mex.cpp meshlib/polyhedron2BSP.cpp meshlib/CPoint.cpp meshlib/CVector.cpp meshlib/Plane3D.cpp meshlib/BSPNode.cpp meshlib/MeshIO.cpp meshlib/FileOperation.cpp meshlib/CPolygon.cpp meshlib/predicates.cpp meshlib/CStopWatch.cpp'])
        eval([mexcommand mexoptions ' -I./meshlib GetListOfConnTri2Tri_mex.cpp meshlib/vector.cpp'])
        eval([mexcommand mexoptions ' -I./meshlib get_ray_shell_intersections.cpp isinvolume_randRay.cpp meshlib/geomath.cpp meshlib/vector.cpp'])
        eval([mexcommand mexoptions ' -I./meshlib expand_bdybuffer_mex.cpp isinvolume_randRay.cpp meshlib/geomath.cpp meshlib/vector.cpp'])
        eval([mexcommand mexoptions ' CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" orient3d_mex.cpp meshlib/geomath.cpp meshlib/vector.cpp meshlib/predicates.cpp -I./meshlib'])
    elseif ~isempty(strfind(os,'GLNX')) % Linux
        eval([mexcommand mexoptions ' -DLINUX -I./meshlib tag_checkerboard3d_mex.cpp meshlib/CStopWatch.cpp'])
        eval([mexcommand mexoptions ' -DLINUX -I./meshlib surface_relations_mex.cpp isinvolume_randRay.cpp meshlib/geomath.cpp meshlib/vector.cpp meshlib/polyhedron2BSP.cpp meshlib/CPoint.cpp meshlib/CVector.cpp meshlib/Plane3D.cpp meshlib/BSPNode.cpp meshlib/MeshIO.cpp meshlib/FileOperation.cpp meshlib/CPolygon.cpp meshlib/predicates.cpp meshlib/CStopWatch.cpp'])
        eval([mexcommand mexoptions ' pnpoly_mex.cpp'])
        eval([mexcommand mexoptions ' -DLINUX -I./meshlib orient_surface_mex.cpp meshlib/vector.cpp meshlib/geomath.cpp isinvolume_randRay.cpp'])
        eval([mexcommand mexoptions ' -DLINUX -I./meshlib involume_mex.cpp isinvolume_randRay.cpp meshlib/geomath.cpp meshlib/vector.cpp'])
        eval([mexcommand mexoptions ' -DLINUX -I./meshlib intersect_ray_shell_mex.cpp ./meshlib/vector.cpp ./meshlib/geomath.cpp'])
        eval([mexcommand mexoptions ' extract_ind_regions_mex.cpp'])
        eval([mexcommand mexoptions ' -DLINUX -I./meshlib PointInPolyhedron_mex.cpp meshlib/polyhedron2BSP.cpp meshlib/CPoint.cpp meshlib/CVector.cpp meshlib/Plane3D.cpp meshlib/BSPNode.cpp meshlib/MeshIO.cpp meshlib/FileOperation.cpp meshlib/CPolygon.cpp meshlib/predicates.cpp meshlib/CStopWatch.cpp'])
        eval([mexcommand mexoptions ' -DLINUX -I./meshlib GetListOfConnTri2Tri_mex.cpp meshlib/vector.cpp'])
        eval([mexcommand mexoptions ' -DLINUX -I./meshlib get_ray_shell_intersections.cpp isinvolume_randRay.cpp meshlib/geomath.cpp meshlib/vector.cpp'])
        eval([mexcommand mexoptions ' -DLINUX -I./meshlib expand_bdybuffer_mex.cpp isinvolume_randRay.cpp meshlib/geomath.cpp meshlib/vector.cpp'])
        eval([mexcommand mexoptions ' -DLINUX CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" orient3d_mex.cpp meshlib/geomath.cpp meshlib/vector.cpp meshlib/predicates.cpp -I./meshlib'])
    end
    
else % PC
    eval([mexcommand mexoptions ' -DWIN32 -I./meshlib tag_checkerboard3d_mex.cpp meshlib/CStopWatch.cpp'])
    eval([mexcommand mexoptions ' -DCPU86 -DWIN32 -I./meshlib surface_relations_mex.cpp isinvolume_randRay.cpp meshlib/geomath.cpp meshlib/vector.cpp meshlib/polyhedron2BSP.cpp meshlib/CPoint.cpp meshlib/CVector.cpp meshlib/Plane3D.cpp meshlib/BSPNode.cpp meshlib/MeshIO.cpp meshlib/FileOperation.cpp meshlib/CPolygon.cpp meshlib/predicates.cpp meshlib/CStopWatch.cpp'])
    eval([mexcommand mexoptions ' pnpoly_mex.cpp'])
    eval([mexcommand mexoptions ' -DWIN32 -DCPU86 -I./meshlib orient_surface_mex.cpp meshlib/vector.cpp meshlib/geomath.cpp isinvolume_randRay.cpp'])
    eval([mexcommand mexoptions ' -DWIN32 -DCPU86 -I./meshlib involume_mex.cpp isinvolume_randRay.cpp meshlib/geomath.cpp meshlib/vector.cpp'])
    eval([mexcommand mexoptions ' -DWIN32 -DCPU86 -I./meshlib intersect_ray_shell_mex.cpp ./meshlib/vector.cpp ./meshlib/geomath.cpp'])
    eval([mexcommand mexoptions ' extract_ind_regions_mex.cpp'])
    eval([mexcommand mexoptions ' -DCPU86 -DWIN32 -I./meshlib PointInPolyhedron_mex.cpp meshlib/polyhedron2BSP.cpp meshlib/CPoint.cpp meshlib/CVector.cpp meshlib/Plane3D.cpp meshlib/BSPNode.cpp meshlib/MeshIO.cpp meshlib/FileOperation.cpp meshlib/CPolygon.cpp meshlib/predicates.cpp meshlib/CStopWatch.cpp'])
    eval([mexcommand mexoptions ' -DWIN32 -DCPU86 -I./meshlib GetListOfConnTri2Tri_mex.cpp meshlib/vector.cpp'])
    eval([mexcommand mexoptions ' -DWIN32 -I./meshlib get_ray_shell_intersections.cpp isinvolume_randRay.cpp meshlib/geomath.cpp meshlib/vector.cpp'])
    eval([mexcommand mexoptions ' -DWIN32 -I./meshlib expand_bdybuffer_mex.cpp isinvolume_randRay.cpp meshlib/geomath.cpp meshlib/vector.cpp'])
    eval([mexcommand mexoptions ' COMPFLAGS="$COMPFLAGS /openmp" LINKFLAGS="$LINKFLAGS /openmp" -DWIN32 -DCPU86 orient3d_mex.cpp meshlib/geomath.cpp meshlib/vector.cpp meshlib/predicates.cpp -I./meshlib'])
end

cd(mypwd)

