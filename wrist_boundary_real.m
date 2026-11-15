function [q, g, h, r] = wrist_boundary_real(p, e, u, time)
%WRIST_BOUNDARY_REAL  Four-electrode boundary conditions for the real
%   PARTICIPANT_1_R wrist cross-section.
%
%   Electrodes lie on the fitted outer skin ellipse:
%     +x arc  : Neumann current injection  n*(sigma*grad V) = g
%     -x arc  : Dirichlet ground            V = 0
%     +y/-y   : insulated (floating measurement electrodes)
%     rest    : insulated
%
%   Total injected current from global WRIST_CURRENT_A (default 1 mA).
%   The current density g is spread uniformly over the electrode arc so the
%   total current equals the requested value exactly.
%
%   Envelope parameters come from global WRIST_GEOM (participant_build_geom).

global WRIST_CURRENT_A WRIST_GEOM
if isempty(WRIST_CURRENT_A)
    WRIST_CURRENT_A = 1e-3;                     % 1 mA
end
if isempty(WRIST_GEOM) || ~isfield(WRIST_GEOM, 'envelope')
    error('wrist_boundary_real:nogeom', 'Global WRIST_GEOM not set.');
end
env = WRIST_GEOM.envelope;
cx = env.cx;  cy = env.cy;
A  = env.aSkin; B = env.bSkin;
halfWidth = 12*pi/180;                          % electrode covers +/- 12 deg

ne = size(e, 2);
q = zeros(1, ne);  g = zeros(1, ne);
h = zeros(1, 2*ne); r = zeros(1, 2*ne);

xm = 0.5*(p(1, e(1, :)) + p(1, e(2, :)));
ym = 0.5*(p(2, e(1, :)) + p(2, e(2, :)));
phi = env.phi;                                  % fitted ellipse rotation
xr =  (xm-cx)*cos(phi) + (ym-cy)*sin(phi);
yr = -(xm-cx)*sin(phi) + (ym-cy)*cos(phi);
rho = sqrt((xr/A).^2 + (yr/B).^2);
ang = atan2(yr, xr);
len = sqrt((p(1, e(1, :)) - p(1, e(2, :))).^2 + (p(2, e(1, :)) - p(2, e(2, :))).^2);

onOuter = rho > 0.995;                          % edges on the skin outer (polygonal) ellipse
winE1 = onOuter & abs(ang)          < halfWidth; % +x inject
winE2 = onOuter & abs(abs(ang)-pi)  < halfWidth; % -x ground

% ---- E1: Neumann injection, total current I spread on the arc ----------
L1 = sum(len(winE1));
if L1 > 0
    g(winE1) = WRIST_CURRENT_A / L1;
end

% ---- E2: Dirichlet ground V = 0 ----------------------------------------
hWin = h(1, 1:ne);
hWin(winE2) = 1;
h(1, 1:ne) = hWin;
h(1, ne+1:2*ne) = hWin;
r(1, 1:2*ne) = 0;
end
