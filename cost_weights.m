function x = cost_weights(p, S, B0, dB, w_prior, w_alpha)
nshapes = size(dB, 1);
T = B0.vertices;
for i=1:nshapes
    T = T + p(i) * dB{i};
end
D = S.vertices - T;
x = sqrt(sum(D.^2, 2));
x = [x; norm(p-w_prior)*w_alpha];
end