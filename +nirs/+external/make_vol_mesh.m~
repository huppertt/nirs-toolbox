clear

load spm_meshes.mat

m{1} = scalp;
m{2} = oskull;
m{3} = iskull;
m{4} = cortex;

bmin = min( scalp.nodes ) - 10;
bmax = max( scalp.nodes ) + 10;

x = bmin(1):2:bmax(1);
y = bmin(2):2:bmax(2);
z = bmin(3):2:bmax(3);

seg = zeros(length(x),length(y),length(z));
for i = 1:length(m)
    oline = logical( surf2vol( m{i}.nodes, m{i}.faces, x, y, z ) );
    
    if i == 4
        sd1 = round(size(oline)/2);
        sd1(1) = round( 0.4*size(oline,1) );
        fill1 = RegionGrowing(oline, 0.5, sd1 );
        
        sd2 = round(size(oline)/2);
        sd2(1) = round( 0.6*size(oline,1) );
        fill2 = RegionGrowing(oline, 0.5, sd2 );
        
        fill = logical(fill1 + fill2);
    else
        fill = RegionGrowing(oline, 0.5, round(size(oline)/2));
    end
    mask{i} = logical( oline+fill );
    
    seg(mask{i}) = i;
end

[n,e,f] = s2m(scalp.nodes, scalp.faces, 0.5, 10);
