function showMeshError(mesh, ref, name)
FV.vertices = mesh.vertices;
FV.faces = mesh.objects.data.vertices;
nverts = size(FV.vertices, 1);
colors = zeros(nverts, 3);
dists = zeros(nverts, 1);
for i=1:nverts
    dists(i) = norm(mesh.vertices(i,:) - ref.vertices(i,:));
end
maxDist = max(dists); minDist = min(dists);
distDiff = maxDist - minDist;
for i=1:nverts
    val = (dists(i)-minDist) / distDiff;
    colors(i, :) = interpolate(1-val);
end

patch(FV, 'FaceVertexCData', colors, 'facecolor', 'interp', 'edgecolor', 'none', 'vertexnormalsmode', 'auto'); camlight;
lighting gouraud;
material dull;
axis vis3d;
axis equal;
caxis([minDist maxDist]);
title(name);
end

function c = interpolate(val)
nrows = 256;
ridx = max(ceil(val*nrows), 1);
cmap = jet(nrows);
c = cmap(ridx,:);
end