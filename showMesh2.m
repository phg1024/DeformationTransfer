function showMesh2(faces, vertices, name)
FV.vertices = vertices;
FV.faces = faces;
patch(FV,'facecolor',[0.5 0.5 0.5], 'edgecolor', 'none', 'vertexnormalsmode', 'auto'); camlight;
lighting gouraud;
material dull;
axis vis3d;
axis equal;
title(name);
end