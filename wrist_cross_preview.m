function wrist_cross_preview(varargin)
%WRIST_CROSS_PREVIEW  Preview PARTICIPANT_1_R wrist cross-sections at
%   several levels so you can pick the z0 for the impedance model.
%
%   wrist_cross_preview
%   wrist_cross_preview('dataDir', '<path to PARTICIPANT_1_R>')
%   wrist_cross_preview('zLevels', [-10 0 10 15], 'outDir', 'results_wrist_real')

cfg = struct( ...
    'dataDir', '', ...
    'zLevels', [-10 -5 0 5 10 15], ...
    'outDir',  'results_wrist_real');
for k = 1:2:numel(varargin)
    if ~isfield(cfg, varargin{k})
        error('wrist_cross_preview:badOption', 'Unknown option: %s', num2str(varargin{k}));
    end
    cfg.(varargin{k}) = varargin{k+1};
end

if ~exist(cfg.outDir, 'dir'), mkdir(cfg.outDir); end
for z0 = cfg.zLevels
    geom = participant_build_geom(cfg.dataDir, z0, 11, true);
    fig = figure('Name', sprintf('z0 = %.1f mm', z0), 'Color', 'w');
    hold on; box on; axis equal;
    th = linspace(0, 2*pi, 300);
    for k = 1:numel(geom.bonePolys)
        fill(geom.bonePolys{k}(:,1), geom.bonePolys{k}(:,2), [0.85 0.65 0.55], 'EdgeColor', 'k');
    end
    for k = 1:numel(geom.cartPolys)
        fill(geom.cartPolys{k}(:,1), geom.cartPolys{k}(:,2), [0.60 0.80 0.85], 'EdgeColor', 'k');
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
        plot(geom.tendons(k,1)+geom.tendons(k,3)*cos(th), geom.tendons(k,2)+geom.tendons(k,3)*sin(th), 'g-');
    end
    if ~isempty(geom.nerve)
        plot(geom.nerve(1)+geom.nerve(3)*cos(th), geom.nerve(2)+geom.nerve(4)*sin(th), 'm-');
    end
    if size(geom.retin,1) >= 3
        fill(geom.retin(:,1), geom.retin(:,2), [0.9 0.9 0.6], 'EdgeColor', 'k');
    end
    xlabel('x [mm]'); ylabel('y [mm]');
    title(sprintf('Wrist cross-section z0 = %.1f mm', z0));
    exportFig(fig, fullfile(cfg.outDir, sprintf('preview_z%+.0f.png', z0)));
    close(fig);
    fprintf('z0 = %+.0f mm : %d bone regions, skin half-axes %.1f x %.1f mm\n', ...
        z0, geom.stats.nBoneRegions, geom.envelope.aSkin, geom.envelope.bSkin);
end
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