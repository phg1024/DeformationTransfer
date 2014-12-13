%% This implementation uses the simplified formulation in Botsch's paper
% Deformation Transfer for Detail-Preserving Surface Editing
function Td = deformationTransfer4(S0, S0grad, T0, T0grad, S1, stationary_indices)
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
%     [v0, v1, v2] = getVertices(S0, S0.faces(i, :));
%     v3 = cross(v1-v0, v2-v0);
%     v3 = v0+v3/norm(v3);
%     V = [v1-v0, v2-v0, v3-v0];
    V = S0grad{i};
    Ds((i-1)*3+1:i*3, 1) = 0.5*dot(cross(V(:,1), V(:,2)), V(:,3));

    [v0t, v1t, v2t] = getVertices(S1, S1.faces(i, :));
    v3t = cross(v1t-v0t, v2t-v0t);
    v3t = v0t+v3t/norm(v3t);
    Vt = [v1t-v0t, v2t-v0t, v3t-v0t];
    S{i} = (Vt/V)';
    
%     [v0, v1, v2] = getVertices(T0, T0.faces(i, :));
%     v3 = cross(v1-v0, v2-v0);
%     v3 = v0+v3/norm(v3);
%     V = [v1-v0, v2-v0, v3-v0];
    V = T0grad{i};
    Vinv = inv(V);
    s = sum(Vinv(1:2, :));
    A = [-s', Vinv(1,:)', Vinv(2,:)'];    
    T{i} = A;
end

c = cell2mat(S);

%assemble matrix A
fprintf('assembling matrix A...\n');
tassemble = tic;
if 0
    A = spalloc(nfaces*3, nverts, 27*nfaces);
    for i=1:nfaces
        verts = S0.faces(i,:);
        rowoffset = (i-1)*3+1;
        Af = T{i};
        for j=1:3
            A(rowoffset:rowoffset+2, verts(j)) = Af(:,j);
        end
    end
else
%     Gi = reshape(repmat([1:nfaces*3], 3, 1), 1, nfaces*9);
%     Gj = reshape(repmat(B{i+1}.faces', 3, 1), 1, nfaces*9);
%     Gs = reshape(cell2mat(T)', 1, nfaces*9);
    Ai = reshape(repmat([1:nfaces*3], 3, 1), 1, nfaces*9);
    Aj = reshape(repmat(S0.faces', 3, 1), 1, nfaces*9);
    As = reshape(cell2mat(T)', 1, nfaces*9);
    A = sparse(Ai, Aj, As, nfaces*3, nverts, 27*nfaces);
end
fprintf('done.\n');
ttotal = toc(tassemble);
fprintf('matrix A assembled in %f seconds\n', ttotal);

% now solve for \tilde{v} = (\tilde v_0, \tilde v_1, \tilde v_2, \tilde v_3)
if 0
    fprintf('post processing matrix A... ');
    A = sparse(A);
    AtD = bsxfun(@times, A', Ds');
    fprintf('done.\n');
    fprintf('solving least square problem ... ');
    tic; x = (AtD*A)\(AtD*c); tsolve = toc;
    fprintf('finished in %f seconds.\n', tsolve);
    Td.vertices = x;
else
    PI = T0.vertices(stationary_indices,:);
    ns = length(stationary_indices);
    II = zeros(ns, nverts);
    II(sub2ind(size(II), 1:ns, stationary_indices')) = 1;
    if 0
        G = sparse(A);
        x = lsqlin(G, c(:,1), [], [], II, PI(:,1));
        y = lsqlin(G, c(:,2), [], [], II, PI(:,2));
        z = lsqlin(G, c(:,3), [], [], II, PI(:,3));
        Td.vertices = [x, y, z];
    else
        G = sparse(A);
        Cmat = [c;PI];
        GtD = bsxfun(@times, G', Ds');        
        G = sparse([G;II]);
        Gt = [GtD, II'];
        Td.vertices = (Gt*G)\(Gt*Cmat);
    end
end
end

function [v0, v1, v2] = getVertices(mesh, face)
v0 = mesh.vertices(face(1),:)';
v1 = mesh.vertices(face(2),:)';
v2 = mesh.vertices(face(3),:)';
end