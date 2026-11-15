function jianmo_real
%JIANMO_REAL  Open PDE Modeler with the PARTICIPANT_1_R wrist cross-section.
%
%   jianmo_real
%   jianmo_real('z0', -10)
%   jianmo_real('dataDir', '<path to PARTICIPANT_1_R>')
%
%   Draws the real bone/cartilage cross-section polygons, the fitted skin/fat
%   envelope and (optional) flexor tendons / median nerve / flexor
%   retinaculum, sets the position-based conductivity coefficient
%   (wrist_sigma_coef_real), and meshes the domain.
%
%   Boundary conditions are NOT auto-applied (edge numbers depend on the
%   geometry decomposition). In the PDE Modeler GUI set the four electrodes
%   manually like before:
%     +x arc  : Neumann,  g = I / L_arc   (e.g. g = 0.5, q = 0)
%     -x arc  : Dirichlet, h = 1, r = 0
%     +/-y    : default insulated (floating measurement)
%   Frequency:  global WRIST_FREQ_HZ;  WRIST_FREQ_HZ = 1e5;

global WRIST_GEOM WRIST_FREQ_HZ
if isempty(WRIST_GEOM)
    WRIST_GEOM = participant_build_geom();
end
if isempty(WRIST_FREQ_HZ)
    WRIST_FREQ_HZ = 1e5;
end
g = WRIST_GEOM;

[pde_fig, ax] = pdeinit;
pdetool('appl_cb', 8);
env = g.envelope;
set(ax, 'DataAspectRatio', [1 1 16]);
set(ax, 'PlotBoxAspectRatio', [1.5 1 1]);
set(ax, 'XLim', [env.cx-env.aSkin-5, env.cx+env.aSkin+5]);
set(ax, 'YLim', [env.cy-env.bSkin-5, env.cy+env.bSkin+5]);
set(ax, 'XTickMode', 'auto');
set(ax, 'YTickMode', 'auto');
pdetool('gridon', 'on');

% ------------------------------------------------------ draw geometry -----
ngd = size(g.gd, 2);
ns = g.ns;
if ischar(ns)                       % classic char matrix: one name per row
    nameOf = @(k) strtrim(ns(k, :));
else                                % cellstr / string vector (R2025b decsg)
    nameOf = @(k) char(ns{k});
end
for k = 1:ngd
    name = nameOf(k);
    col = g.gd(:, k);
    switch col(1)
        case 1                                     % circle
            pdecirc(col(2), col(3), col(4), name);
        case 2                                     % polygon
            n = col(2);
            x = col(3:2+n);
            y = col(3+n:2+2*n);
            pdepoly(x.', y.', name);
        case 4                                     % ellipse
            pdeellip(col(2), col(3), col(4), col(5), col(6), name);
        otherwise
            warning('jianmo_real:type', 'Skipping unsupported gd type %g', col(1));
    end
end
set(findobj(get(pde_fig, 'Children'), 'Tag', 'PDEEval'), 'String', g.sf);

% -------------------------------------------------- PDE coefficients ------
pdeseteq(1, ...
    'wrist_sigma_coef_real', ...
    '0.0', ...
    '0.0', ...
    '1.0', ...
    '0:10', ...
    '0.0', ...
    '0.0', ...
    '[0 100]')
setappdata(pde_fig, 'currparam', char('wrist_sigma_coef_real', '0.0'))

% ---------------------------------------------------------- mesh --------
setappdata(pde_fig, 'trisize', 0.8);
setappdata(pde_fig, 'Hgrad', 1.3);
setappdata(pde_fig, 'refinemethod', 'regular');
setappdata(pde_fig, 'jiggle', char('on', 'mean', ''));
setappdata(pde_fig, 'MesherVersion', 'preR2013a');
pdetool('initmesh')

fprintf('Geometry loaded: real wrist cross-section z0 = %.1f mm\n', g.z0);
fprintf('Next: set 4-electrode boundary conditions manually in the GUI, then Solve.\n');
fprintf('Frequency is controlled by:  global WRIST_FREQ_HZ;  WRIST_FREQ_HZ = 1e5;\n');
end
