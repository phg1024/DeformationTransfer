% function x = cost_weights(p, S, B0, dB)
% nshapes = size(dB, 1);
% T = B0.vertices;
% for i=1:nshapes
%     T = T + p(i) * dB{i};
% end
% D = S.vertices - T;
% x = sqrt(sum(D.^2, 2));
% end

function J = jac_weights(p, S, B0, dB)
nshapes = size(dB, 1);
T = B0.vertices;
for i=1:nshapes
    T = T + p(i) * dB{i};
end
D = S.vertices - T;

nverts = size(T, 1);
J = zeros(nverts, nshapes);
for i=1:nshapes
    J(:,i) = -2.0*sum(D.*dB{i}, 2);
end
end