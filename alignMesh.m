function T0 = alignMesh(T0, T1)

[s, R, t] = estimateTransform(T0.vertices, T1.vertices);
nverts = size(T0.vertices, 1);
for i=1:nverts
    T0.vertices(i,:) = (R * T0.vertices(i,:)' + t)';
end

end