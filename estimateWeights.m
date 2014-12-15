function w = estimateWeights(S, B0, dB, w0, w_prior, w_alpha, maxiters, vis)
if nargin < 7
    maxiters = 128;
end
if nargin < 8
    vis = false;
end

addpath levmar-2.6\matlab\;
nshapes = size(dB, 1);

% initial guess
if nargin < 4
    w = zeros(1, nshapes);
else
    w = w0;
end

nverts = size(B0.vertices, 1);

% measurement vector
meas = zeros(1, nverts+1);

% solve for optimal weights
options=[1e-3, 1E-12, 1E-12, 1E-12, 1E-06];
for i=1:maxiters
    %[ret, popt, info] = levmar('cost_weights', 'jac_weights', w, meas, 5, options, 'unc', S, B0, dB);
    [ret, popt, info] = levmar('cost_weights', 'jac_weights', w, meas, 6, options, 'bc', zeros(1, nshapes), ones(1, nshapes), S, B0, dB, w_prior, w_alpha);
w = popt;
end
fprintf('finished in %d iterations ...\n', ret);

% visualize the fitted mesh
if vis
    T = B0;
    for i=1:nshapes
        T.vertices = T.vertices + popt(i) * dB{i};
    end
    
    figure;showMeshOverlay(T, S, 'overlay');
    figure;showMeshError(T, S, 'error');    
end
end