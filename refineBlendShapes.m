function B_new = refineBlendShapes(S, Sgrad, B, alpha, beta, gamma, prior, w_prior, stationary_indices)

nverts = size(S{1}.vertices, 1);
nfaces = size(S{1}.faces, 1);
nposes = size(S, 1);
nshapes = size(B, 1)-1;

B0 = B{1};

Bmat = zeros(1, 9*nfaces);
Btri = cell(nfaces, 1);
% triangle areas of B0
D = zeros(nfaces*3, 1);
% assemble the matrices
tic;
for j=1:nfaces
    jstart = (j-1)*9+1; jend = j*9;
    [Bj, dj] = triangleGradient(B0, j);
    Bmat(1, jstart:jend) = reshape(Bj, 1, 9);
    Btri{j} = reshape(Bj, 1, 9);
    D((j-1)*3+1:j*3, 1) = dj;
end

A = cell2mat(alpha);
toc

tic;
Mtri = cell(1, nfaces);
% face by face
parfor j=1:nfaces
    % all poses, only one face
    Smat = Sgrad{j};
    jstart = (j-1)*9+1; jend = j*9;
    Rhs_tri = Smat - repmat(Btri{j}, nposes, 1);
    % augment A and Rhs with prior
    A_prior = beta*diag(w_prior(:,j));
    Amat = [A; A_prior];
    Rhs_prior = A_prior * prior(:,jstart:jend);
    Rhs = [Rhs_tri; Rhs_prior];
    % Mmat is the deformation gradients of all the blendshapes
    Mtri{j} = (Amat'*Amat)\(Amat'*Rhs);
end
Mmat = cell2mat(Mtri);
toc

% reconstruct the blendshapes
B_new = cell(nshapes,1);
parfor i=1:nshapes
    fprintf('reconstructing blendshape %d...\n', i);
    M0 = reshape(Bmat, 9, nfaces);
    Mi = reshape(Mmat(i,:), 9, nfaces);
    Cmat = zeros(3*nfaces, 3);
    T = cell(nfaces, 1);    
    for j=1:nfaces
        M0j = reshape(M0(:,j), 3, 3);
        Mij = reshape(Mi(:,j), 3, 3) + M0j;
        
        r0 = (j-1)*3+1; r1 = j*3;        
        Cmat(r0:r1,:) = (Mij/M0j)';
        
        M0j_inv = inv(M0j);
        s = sum(M0j_inv(1:2, :));
        T{j} = [-s', M0j_inv(1,:)', M0j_inv(2,:)'];  
    end
    
    if 1
        Gi = reshape(repmat([1:nfaces*3], 3, 1), 1, nfaces*9);
        Gj = reshape(repmat(B{i+1}.faces', 3, 1), 1, nfaces*9);
        Gs = reshape(cell2mat(T)', 1, nfaces*9);
        G = sparse(Gi, Gj, Gs, nfaces*3, nverts, 27*nfaces);
    else
        % Gi    Gj      Gs
        % 1     vj1     Tj_11
        % 1     vj2     Tj_12
        % 1     vj3     Tj_13
        % 2     vj1
        % 2     vj2
        % 2     vj3
        % 3     vj1
        % 3     vj2
        % 3     vj3
        G = spalloc(nfaces*3, nverts, 27*nfaces);
        for j=1:nfaces
            verts = B{i+1}.faces(j,:);
            rowoffset = (j-1)*3+1;
            Af = T{j};
            for k=1:3
                G(rowoffset:rowoffset+2, verts(k)) = Af(:,k);
            end
        end
    end

    % now solve for \tilde{v} = (\tilde v_0, \tilde v_1, \tilde v_2, \tilde v_3)
    G = sparse(G);
    %GtD = bsxfun(@times, G', D');
    %x = (GtD*G)\(GtD*Cmat);
    PI = B0.vertices(stationary_indices,:);
    ns = length(stationary_indices);
    II = zeros(ns, nverts);
    II(sub2ind(size(II), 1:ns, stationary_indices')) = 1;
    if 0
        x = lsqlin(G, Cmat(:,1), [], [], II, PI(:,1));
        y = lsqlin(G, Cmat(:,2), [], [], II, PI(:,2));
        z = lsqlin(G, Cmat(:,3), [], [], II, PI(:,3));
        B_new{i} = [x, y, z];
    else
        Mrhs = [Cmat;PI];
        GtD = sparse(bsxfun(@times, G', D'));
        G = sparse([G;II]);
        Gt = sparse([GtD, sparse(II'*gamma)]);
        B_new{i} = (Gt*G)\(Gt*Mrhs);
    end
end

end