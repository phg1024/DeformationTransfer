function showMeshError(mesh, ref, name)
FV.vertices = mesh.vertices;
FV.faces = mesh.objects.data.vertices;
nverts = size(FV.vertices, 1);
colors = zeros(nverts, 3);
dists = zeros(nverts, 1);

aligned = 1;
if aligned
    center1 = zeros(1, 3);
    center2 = zeros(1, 3);
else
    center1 = mean(mesh.vertices);
    center2 = mean(ref.vertices);
end

for i=1:nverts
    dists(i) = norm((mesh.vertices(i,:)-center1) - (ref.vertices(i,:)-center2));
end
maxDist = max(dists); minDist = min(dists);
distDiff = maxDist - minDist;
for i=1:nverts
    val = (dists(i)-minDist) / distDiff;
    colors(i, :) = interpolate(val);
end

patch(FV, 'FaceVertexCData', colors, 'facecolor', 'interp', 'edgecolor', 'none', ...
    'vertexnormalsmode', 'auto'); 
camlight;
lighting gouraud;
material dull;
axis vis3d;
axis equal;
colorbar;
colormap jet(256);
caxis([minDist maxDist]);
title(name);
end

function c = interpolate(val)
nrows = 256;
ridx = max(ceil(val*nrows), 1);
cmap = jet(nrows);
c = cmap(ridx,:);
end