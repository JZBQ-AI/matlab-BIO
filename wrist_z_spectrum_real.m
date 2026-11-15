function wrist_z_spectrum_real(varargin)
%WRIST_Z_SPECTRUM_REAL  Frequency sweep on the PARTICIPANT_1_R wrist
%   cross-section model (2-D PDE FEM, same solver family as wrist_z_spectrum).
%
%   wrist_z_spectrum_real
%   wrist_z_spectrum_real('z0', -10)                % distal forearm level
%   wrist_z_spectrum_real('z0', 10, 'Hmax', 0.6)
%   wrist_z_spectrum_real('freqs', logspace(2, 5, 13))
%   wrist_z_spectrum_real('complex', false)         % sigma only (real solve)
%
%   Geometry: real bone/cartilage cross-sections from PARTICIPANT_1_R at the
%   chosen z0, plus a fitted skin/fat envelope and (optional) approximate
%   flexor tendons, median nerve and flexor retinaculum.
%
%   Electrodes (same convention as before):
%     +x arc : inject I,  -x arc : ground,  +/-y arcs : measurement (floating)
%
%   Outputs: console table, CSV, four-panel spectrum figure, cross-section
%   figure, field figure, and a .mat with geom + gd/ns/sf for PDE Modeler.
%
%   Requires: Partial Differential Equation Toolbox (decsg/initmesh/assema/
%   assempde). Files: participant_read_inp, participant_cross_section,
%   participant_build_geom, wrist_tissue_sigma_real, wrist_boundary_real.

cfg = struct( ...
    'dataDir',    '', ...
    'z0',         10, ...
    'marginMm',   11, ...
    'includeSoft', true, ...
    'freqs',      logspace(2, 6, 25), ...
    'Hmax',       0.8, ...
    'Hgrad',      1.3, ...
    'I',          1e-3, ...
    'complex',    true, ...
    'outDir',     'results_wrist_real', ...
    'saveCsv',    'wrist_z_spectrum_real.csv', ...
    'makePlot',   true);
for k = 1:2:numel(varargin)
    if ~isfield(cfg, varargin{k})
        error('wrist_z_spectrum_real:badOption', 'Unknown option: %s', num2str(varargin{k}));
    end
    cfg.(varargin{k}) = varargin{k+1};
end
cfg.freqs = cfg.freqs(:).';

% ============================================================ 1) geometry ====
fprintf('\n================ PARTICIPANT_1_R wrist model ================\n');
geom = participant_build_geom(cfg.dataDir, cfg.z0, cfg.marginMm, cfg.includeSoft);
fprintf('z0 = %.1f mm | skin half-axes %.1f x %.1f mm | bone span %.1f x %.1f mm\n', ...
    geom.z0, geom.envelope.aSkin, geom.envelope.bSkin, ...
    geom.stats.boneSpanX_mm, geom.stats.boneSpanY_mm);
fprintf('bone regions: %d, cartilage regions: %d\n', ...
    geom.stats.nBoneRegions, geom.stats.nCartRegions);

global WRIST_GEOM WRIST_CURRENT_A
WRIST_GEOM = geom;
WRIST_CURRENT_A = cfg.I;      % used by wrist_boundary_real for the BC matrices

% ------------------------------------------------------------ mesh ---------
dl = decsg(geom.gd, geom.sf, geom.ns);
[p, e, t] = meshDomain(dl, cfg.Hmax, cfg.Hgrad);
np = size(p, 2);
fprintf('mesh: %d nodes, %d elements\n', np, size(t, 2));

% ------------------------------------------------------------ boundary -----
[~, ~, ~, Q, G, H, R] = assempde(@wrist_boundary_real, p, e, t, 1, 0, 0);

% electrode nodes on the outer skin ellipse
env = geom.envelope;
phi = env.phi;                                  % fitted ellipse rotation
xrP =  (p(1,:)-env.cx)*cos(phi) + (p(2,:)-env.cy)*sin(phi);
yrP = -(p(1,:)-env.cx)*sin(phi) + (p(2,:)-env.cy)*cos(phi);
rhoP = sqrt((xrP/env.aSkin).^2 + (yrP/env.bSkin).^2);
angP = atan2(yrP, xrP);
onOut = rhoP > 0.995;   % outer skin is a 128-gon, allow small chord inset
halfW = 12*pi/180;
nodeE1 = onOut & abs(angP)         < halfW;   % +x inject
nodeE2 = onOut & abs(abs(angP)-pi) < halfW;   % -x ground
nodeE3 = onOut & abs(angP-pi/2)    < halfW;   % +y measure
nodeE4 = onOut & abs(angP+pi/2)    < halfW;   % -y measure

% element centroids -> tissue assignment
xc = mean(reshape(p(1, t(1:3, :)), 3, []), 1);
yc = mean(reshape(p(2, t(1:3, :)), 3, []), 1);
nt = numel(xc);

% ============================================================ 2) sweep ======
freqs = cfg.freqs;
nf = numel(freqs);
Zt = zeros(1, nf);
Zb = zeros(1, nf);
complexOK = cfg.complex;

fprintf('\n%-12s %12s %12s %12s %12s\n', 'f [Hz]', '|Zt| [ohm]', 'phase [deg]', 'R [ohm]', 'X [ohm]');
for k = 1:nf
    f = freqs(k);
    [sigTri, eprTri] = wrist_tissue_sigma_real(xc, yc, f);   % 1 x nt

    if complexOK
        cTri = sigTri + 1i*2*pi*f*8.8541878128e-12*eprTri;
    else
        cTri = sigTri;
    end

    K = sparse(np, np);
    uv = unique(cTri);
    M0 = sparse(np, np);       % empty mass matrix (see note below)
    try
        for m = 1:numel(uv)
            idx = find(cTri == uv(m));
            if isempty(idx), continue; end
            [ks, ~, ~] = assema(p, t(:, idx), uv(m), 0, 0);
            K = K + ks;
        end
        % NB: the 2nd argument (M) must be an empty *sparse* matrix. Passing
        % the scalar 0 makes assempde evaluate K+M+Q, and "sparse + scalar"
        % densifies the 27k x 27k matrix (~150 s and several GB per solve).
        u = assempde(K, M0, zeros(np, 1), Q, G, H, R);
    catch
        if complexOK
            warning('wrist_z_spectrum_real:complex', ...
                ['Complex coefficients are not supported by this toolbox ' ...
                 'version. Falling back to real (sigma-only) solve.']);
            complexOK = false;
            cTri = sigTri;
            uv2 = unique(cTri);
            K = sparse(np, np);
            for m = 1:numel(uv2)
                idx = find(cTri == uv2(m));
                if isempty(idx), continue; end
                [ks, ~, ~] = assema(p, t(:, idx), uv2(m), 0, 0);
                K = K + ks;
            end
            u = assempde(K, M0, zeros(np, 1), Q, G, H, R);
        else
            rethrow(lasterror);
        end
    end

    V1 = mean(u(nodeE1));  V2 = mean(u(nodeE2));
    V3 = mean(u(nodeE3));  V4 = mean(u(nodeE4));
    Zt(k) = (V3 - V4) / cfg.I;
    Zb(k) = (V1 - V2) / cfg.I;

    fprintf('%-12.4g %12.4f %12.4f %12.4f %12.4f\n', ...
        f, abs(Zt(k)), angle(Zt(k))*180/pi, real(Zt(k)), imag(Zt(k)));

    if k == ceil(nf/2)
        Vfield = u; fField = f;   % keep one field solution for plotting
    end
end

% ============================================================ 3) output =====
if ~exist(cfg.outDir, 'dir'), mkdir(cfg.outDir); end

T = table(freqs(:), abs(Zt).', angle(Zt).'*180/pi, real(Zt).', imag(Zt).', abs(Zb).', ...
    'VariableNames', {'f_Hz','Zt_mag_ohm','Zt_phase_deg','Zt_real_ohm','Zt_imag_ohm','Zb_mag_ohm'});
writetable(T, fullfile(cfg.outDir, cfg.saveCsv));
gd = geom.gd;      %#ok<NASGU>
ns = geom.ns;      %#ok<NASGU>
sf = geom.sf;      %#ok<NASGU>
save(fullfile(cfg.outDir, 'wrist_geom_real.mat'), 'geom', 'gd', 'ns', 'sf', 'p', 'e', 't', 'T', 'Zt', 'Zb', 'freqs');
fprintf('\nresults saved: %s\n', fullfile(cfg.outDir, cfg.saveCsv));

if cfg.makePlot
    plotCrossSection(geom, cfg);
    plotSpectrum(freqs, Zt, Zb, cfg);
    if exist('Vfield', 'var')
        plotField(p, e, t, Vfield, fField, cfg);
    end
end
fprintf('==============================================================\n');
end


% ==================================================================== plots ==
function plotCrossSection(geom, cfg)
fig = figure('Name', sprintf('Real wrist cross-section z0 = %.1f mm', geom.z0), 'Color', 'w');
hold on; box on; axis equal;
th = linspace(0, 2*pi, 300);
for k = 1:numel(geom.bonePolys)
    fill(geom.bonePolys{k}(:,1), geom.bonePolys{k}(:,2), [0.85 0.65 0.55], 'EdgeColor', 'k');
end
for k = 1:numel(geom.cartPolys)
    fill(geom.cartPolys{k}(:,1), geom.cartPolys{k}(:,2), [0.6 0.8 0.85], 'EdgeColor', 'k');
end
env = geom.envelope;
cp = cos(env.phi); sp = sin(env.phi);
plot(env.cx + env.aSkin*cos(th)*cp - env.bSkin*sin(th)*sp, ...
     env.cy + env.aSkin*cos(th)*sp + env.bSkin*sin(th)*cp, 'k-', 'LineWidth', 1.5);
plot(env.cx + env.aSkinIn*cos(th)*cp - env.bSkinIn*sin(th)*sp, ...
     env.cy + env.aSkinIn*cos(th)*sp + env.bSkinIn*sin(th)*cp, 'k--');
plot(env.cx + env.aFatIn*cos(th)*cp - env.bFatIn*sin(th)*sp, ...
     env.cy + env.aFatIn*cos(th)*sp + env.bFatIn*sin(th)*cp, 'k--');
for k = 1:size(geom.tendons, 1)
    plot(geom.tendons(k,1) + geom.tendons(k,3)*cos(th), geom.tendons(k,2) + geom.tendons(k,3)*sin(th), 'g-');
end
if ~isempty(geom.nerve)
    plot(geom.nerve(1) + geom.nerve(3)*cos(th), geom.nerve(2) + geom.nerve(4)*sin(th), 'm-');
end
if size(geom.retin, 1) >= 3
    fill(geom.retin(:,1), geom.retin(:,2), [0.9 0.9 0.6], 'EdgeColor', 'k');
end
xlabel('x [mm]'); ylabel('y [mm]');
title(sprintf('Real wrist cross-section z0 = %.1f mm (bone=brown, cartilage=cyan, tendon=green, nerve=magenta, retinaculum=yellow)', geom.z0));
if ~exist(cfg.outDir, 'dir'), mkdir(cfg.outDir); end
exportFig(fig, fullfile(cfg.outDir, 'cross_section.png'));
end


function plotSpectrum(freqs, Zt, Zb, cfg)
fig = figure('Name', 'Real-geometry wrist impedance spectrum', 'Color', 'w');
subplot(2,2,1); semilogx(freqs, abs(Zt), 'o-', 'LineWidth', 1.5); grid on
xlabel('f [Hz]'); ylabel('|Z_t| [\Omega]'); title('4-electrode |Z|');
subplot(2,2,2); semilogx(freqs, angle(Zt)*180/pi, 'o-', 'LineWidth', 1.5); grid on
xlabel('f [Hz]'); ylabel('phase [deg]'); title('4-electrode phase');
subplot(2,2,3); plot(real(Zt), -imag(Zt), 'o-', 'LineWidth', 1.5); grid on; axis equal
xlabel('R [\Omega]'); ylabel('-X [\Omega]'); title('Nyquist');
subplot(2,2,4); semilogx(freqs, abs(Zb), 's-', 'LineWidth', 1.5); grid on
xlabel('f [Hz]'); ylabel('|Z_b| [\Omega]'); title('2-electrode |Z|');
if ~exist(cfg.outDir, 'dir'), mkdir(cfg.outDir); end
exportFig(fig, fullfile(cfg.outDir, 'impedance_spectrum.png'));
end


function plotField(p, e, t, u, f, cfg)
fig = figure('Name', sprintf('Potential field f = %.0f Hz', f), 'Color', 'w');
try
    pdeplot(p, e, t, 'XYData', real(u), 'ZData', real(u));
    colorbar;
catch
    trisurf(t(1:3,:).', p(1,:), p(2,:), real(u));
    colorbar;
end
xlabel('x [mm]'); ylabel('y [mm]'); title(sprintf('Re(V) [V], f = %.0f Hz', f));
if ~exist(cfg.outDir, 'dir'), mkdir(cfg.outDir); end
exportFig(fig, fullfile(cfg.outDir, sprintf('field_%gHz.png', f)));
end
function exportFig(fig, path)
%EXPORTFIG  Save figure as PNG on both new releases (exportgraphics) and
%   older releases that only provide print.
if exist('exportgraphics', 'file')
    exportgraphics(fig, path, 'Resolution', 150);
else
    folder = fileparts(path);
    if ~isempty(folder) && ~exist(folder, 'dir')
        mkdir(folder);
    end
    print(fig, path, '-dpng', '-r150');
end
end


function [p, e, t] = meshDomain(dl, Hmax, Hgrad)
%MESHDOMAIN  initmesh with a fallback to the pre-R2013a mesher.
%   The R2025b mesher is stricter about short edges; if it refuses a layer
%   (e.g. a cartilage sliver), the classic mesher still meshes it.
try
    [p, e, t] = initmesh(dl, 'Hmax', Hmax, 'Hgrad', Hgrad);
catch ME
    warning('wrist_z_spectrum_real:mesher', ...
        'Default mesher failed (%s); retrying with MesherVersion R2013a.', ME.message);
    [p, e, t] = initmesh(dl, 'Hmax', Hmax, 'Hgrad', Hgrad, 'MesherVersion', 'R2013a');
end
end