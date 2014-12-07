function showMesh(mesh, name)
FV.vertices = mesh.vertices;
FV.faces = mesh.faces;
patch(FV,'facecolor',[0.5 0.5 0.5], 'edgecolor', 'none', 'vertexnormalsmode', 'auto'); camlight;
lighting gouraud;
material dull;
axis vis3d;
axis equal;
title(name);
end