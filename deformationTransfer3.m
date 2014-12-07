%% This implementation uses the simplified formulation in Botsch's paper
% Deformation Transfer for Detail-Preserving Surface Editing
function Td = deformationTransfer3(S0, T0, S1)
% make a copy of T0
Td = T0;

% solve it with vertex formulation

% for each face in the source mesh, compute the deformation matrix
nfaces = size(S0.faces, 1);
nverts = size(S0.vertices, 1);

S = cell(nfaces, 1);
T = cell(nfaces, 1);
Ds = zeros(3*nfaces, 1);   % triangle areas
for i=1:nfaces
    [v0, v1, v2] = getVertices(S0, S0.faces(i, :));
    v3 = cross(v1-v0, v2-v0);
    v3 = v0+v3/norm(v3);
    V = [v1-v0, v2-v0, v3-v0];

    [v0t, v1t, v2t] = getVertices(S1, S1.faces(i, :));
    v3t = cross(v1t-v0t, v2t-v0t);
    v3t = v0t+v3t/norm(v3t);
    Vt = [v1t-v0t, v2t-v0t, v3t-v0t];
    S{i} = (Vt/V)';
    Ds((i-1)*3+1:i*3, 1) = 0.5*dot(cross(v1-v0, v2-v0), v3-v0);
    
    [v0, v1, v2] = getVertices(T0, T0.faces(i, :));
    v3 = cross(v1-v0, v2-v0);
    v3 = v0+v3/norm(v3);
    V = [v1-v0, v2-v0, v3-v0];
    Vinv = inv(V);
    s = sum(Vinv(1:2, :));
    A = [-s', Vinv(1,:)', Vinv(2,:)'];    
    T{i} = A;
end

c = cell2mat(S);

%assemble matrix A
fprintf('assembling matrix A...\n');
A = spalloc(nfaces*3, nverts, 27*nfaces);
ttotal = 0;
tic; tic;
for i=1:nfaces
    if mod(i,10000) == 0
        tstep = toc;
        ttotal = ttotal + tstep;
        fprintf('processed %d faces in %f seconds\n', i, ttotal);        
        tic;
    end
    verts = S0.faces(i,:);
    rowoffset = (i-1)*3+1;
    Af = T{i};
    for j=1:3
        A(rowoffset:rowoffset+2, verts(j)) = Af(:,j);
    end
end
toc;
fprintf('done.\n');
ttotal = ttotal + toc;
fprintf('matrix A assembled in %f seconds\n', ttotal);

% now solve for \tilde{v} = (\tilde v_0, \tilde v_1, \tilde v_2, \tilde v_3)
fprintf('post processing matrix A...\n');
A = sparse(A);
AtD = bsxfun(@times, A', Ds');
fprintf('done.\n');
fprintf('solving least square problem ...\n');
tic; x = (AtD*A)\(AtD*c); tsolve = toc;
fprintf('solved in %f seconds.\n', tsolve);
Td.vertices = x;
end

function [v0, v1, v2] = getVertices(mesh, face)
v0 = mesh.vertices(face(1),:)';
v1 = mesh.vertices(face(2),:)';
v2 = mesh.vertices(face(3),:)';
end