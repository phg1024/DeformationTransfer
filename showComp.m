figure;
subplot(2, 3, 1), showMesh(S0, 'source');
subplot(2, 3, 2), showMesh(S1, 'deformed');
subplot(2, 3, 4), showMesh(T0, 'target');
subplot(2, 3, 5), showMesh(T1_tf, 'transferred');
subplot(2, 3, 6), showMesh(T1, 'reference');
subplot(2, 3, 3), showMeshError(T1, T1_tf, 'error'); colorbar;
figure;
showMeshError(T1, T1_tf, 'error'); colorbar;