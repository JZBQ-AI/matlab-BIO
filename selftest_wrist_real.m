function selftest_wrist_real(dataDir)
%SELFTEST_WRIST_REAL  Quick end-to-end check of the real wrist model.
%   selftest_wrist_real
%   selftest_wrist_real('<full path to PARTICIPANT_1_R>')
%
%   Checks: cross-sections at z0 = 10 and -10 mm, the zero-padded polygon gd
%   matrix (all columns simple, no repeated vertices), decsg/initmesh (with the
%   R2013a mesher as a fallback), and one 100 kHz PDE solve.

if nargin < 1 || isempty(dataDir)
    dataDir = '';
end

fprintf('--- [1/4] build geometry at z0 = +10 mm (carpal row) ---\n');
g10 = participant_build_geom(dataDir, 10);
fprintf('bone regions = %d, cartilage regions = %d, gd size = %d x %d\n', ...
    numel(g10.bonePolys), numel(g10.cartPolys), size(g10.gd, 1), size(g10.gd, 2));
for k = 1:numel(g10.bonePolys)
    a = abs(polyarea(g10.bonePolys{k}(:, 1), g10.bonePolys{k}(:, 2)));
    fprintf('  %-30s area = %8.1f mm^2, points = %d\n', ...
        g10.boneNames{k}, a, size(g10.bonePolys{k}, 1));
end
fprintf('skin envelope: center (%.2f, %.2f), half-axes %.1f x %.1f mm\n', ...
    g10.envelope.cx, g10.envelope.cy, g10.envelope.aSkin, g10.envelope.bSkin);

fprintf('--- [2/4] decsg + mesh ---\n');
dl = decsg(g10.gd, g10.sf, g10.ns);
fprintf('geometric repairs: %d ring(s) rebuilt after the plane cut\n', numel(g10.repairLog));
[p, e, t] = meshDomain(dl, 1.2, 1.3);
np = size(p, 2);
fprintf('mesh: %d nodes, %d triangles\n', np, size(t, 2));

fprintf('--- [3/4] single PDE solve at 100 kHz (real conductivity) ---\n');
global WRIST_GEOM WRIST_CURRENT_A
WRIST_GEOM = g10;
WRIST_CURRENT_A = 1e-3;
xc = mean(reshape(p(1, t(1:3, :)), 3, []), 1);
yc = mean(reshape(p(2, t(1:3, :)), 3, []), 1);
sig = wrist_tissue_sigma_real(xc, yc, 1e5);
K = sparse(np, np);
uv = unique(sig);
for m = 1:numel(uv)
    idx = find(sig == uv(m));
    if isempty(idx)
        continue;
    end
    [ks, ~, ~] = assema(p, t(:, idx), uv(m), 0, 0);
    K = K + ks;
end
[~, ~, ~, Q, G, H, R] = assempde(@wrist_boundary_real, p, e, t, 1, 0, 0);
% NB: pass an empty *sparse* mass matrix, not the scalar 0. assempde computes
% K+M+Q internally and "sparse + scalar" densifies the 27k x 27k matrix
% (that alone took ~150 s and several GB before this was fixed).
u = assempde(K, sparse(np, np), zeros(np, 1), Q, G, H, R);
fprintf('solve OK: potential range [%.4f, %.4f] V for 1 mA injection\n', ...
    min(real(u)), max(real(u)));

fprintf('--- [4/4] build geometry at z0 = -10 mm (distal forearm) ---\n');
gm10 = participant_build_geom(dataDir, -10);
fprintf('bone regions = %d (%s)\n', numel(gm10.bonePolys), strjoin(gm10.boneNames, ', '));

fprintf('\nSELF TEST PASSED. Now run: run_wrist_real\n');
end


function [p, e, t] = meshDomain(dl, Hmax, Hgrad)
%MESHDOMAIN  initmesh with a fallback to the pre-R2013a mesher.
%   The R2025b mesher is stricter about short edges; if it refuses a layer,
%   the classic mesher still produces a usable triangulation.
try
    [p, e, t] = initmesh(dl, 'Hmax', Hmax, 'Hgrad', Hgrad);
catch ME
    warning('selftest_wrist_real:mesher', ...
        'Default mesher failed (%s); retrying with MesherVersion R2013a.', ME.message);
    [p, e, t] = initmesh(dl, 'Hmax', Hmax, 'Hgrad', Hgrad, 'MesherVersion', 'R2013a');
end
end
