%% only works for quad mesh
function M = triangulateMesh(S)
M = S;
faces = S.faces;
nfaces = size(faces, 1);
M.faces = zeros(nfaces*2, 3);
M.faces(1:nfaces,:) = S.faces(:,[1 2 3]);
M.faces(nfaces+1:end,:) = S.faces(:,[1 3 4]);
end