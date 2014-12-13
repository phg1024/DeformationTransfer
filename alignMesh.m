function T0 = alignMesh(T0, T1, indices)

if nargin < 3
    nverts = size(T0.vertices, 1);
    indices = [1:nverts]';
end
[s, R, t] = estimateTransform(T0.vertices(indices,:), T1.vertices(indices,:));
nverts = size(T0.vertices, 1);
% for i=1:nverts
%     T0.vertices(i,:) = (R * T0.vertices(i,:)' + t)';
% end
T0.vertices = T0.vertices + repmat(t', nverts, 1);

end