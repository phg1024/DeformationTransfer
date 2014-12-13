% formulate this as a linear optimization
% [v1 v2 ... vn] * [w1 w2 ... wn]' = vs - v0
function w = estimateWeights_linear(S, B0, dB, stationary_indices, vis)
if nargin < 5
    vis = false;
end

nshapes = size(dB, 1);

nverts = size(B0.vertices, 1);

V = zeros(nverts*3, nshapes);
for j=1:nshapes
    V(:,j) = reshape(dB{j}, nverts*3, 1);
end

Rhs = reshape(S.vertices - B0.vertices, nverts*3, 1);

w = (V'*V)\(V'*Rhs);
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