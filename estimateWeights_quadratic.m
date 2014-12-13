% formulate this as a linear optimization
% [v1 v2 ... vn] * [w1 w2 ... wn]' = vs - v0
function w = estimateWeights_quadratic(S, B0, dB, w0, stationary_indices, vis)
if nargin < 5
    vis = false;
end

nshapes = size(dB, 1);

nverts = size(B0.vertices, 1);

H = zeros(nshapes, nshapes);
f = zeros(1, nshapes);
for j=1:nverts
    V = zeros(3, nshapes);
    for i=1:nshapes
        V(:,i) = dB{i}(j,:)';
    end
    H = H + V'*V;
    
    vsv0 = S.vertices(i,:) - B0.vertices(i,:);
    f = f + vsv0 * V;
end

w = quadprog(H, -2.0*f, [], [], [], [], zeros(1, nshapes), ones(1, nshapes), w0);
w = w';

% visualize the fitted mesh
if vis
    T = B0;
    for i=1:nshapes
        T.vertices = T.vertices + w(i) * dB{i};
    end
    T = alignMesh(T, S, stationary_indices);
    figure;showMeshOverlay(T, S, 'overlay');
    figure;showMeshError(T, S, 'error');    
end
end