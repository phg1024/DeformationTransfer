function Td = deformationTransfer(S0, T0, S1)
% make a copy of T0
Td = T0;

% solve it with vertex formulation

% for each face in the source mesh, compute the deformation matrix
nfaces = size(S0.faces, 1);
nverts = size(S0.vertices, 1);

% the 
S = cell(nfaces, 1);
T = cell(nfaces, 1);
for i=1:nfaces
    [v0, v1, v2] = getVertices(S0, S0.faces(i, :));
    v3 = cross(v1-v0, v2-v0);
    v3 = v0+v3/norm(v3);
    V = [v1-v0, v2-v0, v3-v0];

    [v0t, v1t, v2t] = getVertices(S1, S1.faces(i, :));
    v3t = cross(v1t-v0t, v2t-v0t);
    v3t = v0t+v3t/norm(v3t);
    Vt = [v1t-v0t, v2t-v0t, v3t-v0t];
    S{i} = Vt / V;
    
    [v0, v1, v2] = getVertices(T0, T0.faces(i, :));
    v3 = cross(v1-v0, v2-v0);
    v3 = v0+v3/norm(v3);
    V = [v1-v0, v2-v0, v3-v0];
    Vinv = inv(V);
    s = sum(Vinv);
    A = zeros(9, 12);
    A(1:3,:) = [-s(1) 0 0 Vinv(1, 1) 0 0 Vinv(2, 1) 0 0 Vinv(3, 1) 0 0; ...
         -s(2) 0 0 Vinv(1, 2) 0 0 Vinv(2, 2) 0 0 Vinv(3, 2) 0 0; ...
         -s(3) 0 0 Vinv(1, 3) 0 0 Vinv(2, 3) 0 0 Vinv(3, 3) 0 0;];
    A(4:6,:) = [zeros(3, 1), A(1:3,1:end-1)];
    A(7:9,:) = [zeros(3, 2), A(1:3,1:end-2)];
    T{i} = A;
end

c = reshape(cell2mat(S)', nfaces*9, 1);

%assemble matrix A
fprintf('assembling matrix A...\n');
A = sparse(nfaces*9, (nverts+nfaces)*3);
for i=1:nfaces
    if mod(i,1000) == 0
        fprintf('processed %d faces\n', i);
    end
    verts = S0.faces(i,:);
    rowoffset = (i-1)*9+1;
    Af = T{i};
    for j=1:3
        vidx = verts(j);
        coloffset = (vidx-1) * 3 + 1;
        A(rowoffset:rowoffset+8, coloffset:coloffset+2) = Af(:, (j-1)*3+1:j*3);
    end
    A(rowoffset:rowoffset+8, (nverts+i-1)*3+1:(nverts+i)*3) = Af(:, 10:12);
end
fprintf('done.\n');

% now solve for \tilde{v} = (\tilde v_0, \tilde v_1, \tilde v_2, \tilde v_3)
fprintf('solving least square problem ...\n');
A = sparse(A);
x = (A'*A)\(A'*c);
fprintf('done.\n');
size(x)
x = reshape(x, 3, nverts+nfaces);
Td.vertices = x(:,1:nverts)';
end

function [v0, v1, v2] = getVertices(mesh, face)
v0 = mesh.vertices(face(1),:)';
v1 = mesh.vertices(face(2),:)';
v2 = mesh.vertices(face(3),:)';
end