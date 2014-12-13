function [M, d] = triangleGradient(mesh, fidx)
[v0, v1, v2] = getVertices(mesh, fidx);
n = cross(v1-v0, v2-v0);
n = n/norm(n);
M = [v1-v0, v2-v0, n];
d = 0.5 * dot(cross(v1-v0, v2-v0), n);
end

function [v0, v1, v2] = getVertices(mesh, fidx)
face = mesh.faces(fidx,:);
v0 = mesh.vertices(face(1),:)';
v1 = mesh.vertices(face(2),:)';
v2 = mesh.vertices(face(3),:)';
end