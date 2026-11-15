function [sig, epr] = wrist_tissue_sigma_real(x, y, f)
%WRIST_TISSUE_SIGMA_REAL  Position-based tissue conductivity/permittivity for
%   the PARTICIPANT_1_R wrist cross-section model.
%
%   [sig, epr] = wrist_tissue_sigma_real(x, y, f)
%
%   x, y : coordinates in mm (row vectors, as used by PDE coefficient calls)
%   f    : frequency in Hz
%   sig  : conductivity [S/m]
%   epr  : relative permittivity [-]
%
%   Requires global WRIST_GEOM, a struct from participant_build_geom, set by
%   wrist_z_spectrum_real or jianmo_real before solving.
%
%   Tissue priority: bone > cartilage > nerve > tendon > retinaculum >
%   skin > fat > deep soft tissue.

global WRIST_GEOM
if isempty(WRIST_GEOM) || ~isfield(WRIST_GEOM, 'bonePolys')
    error('wrist_tissue_sigma_real:nogeom', ...
        sprintf(['Global WRIST_GEOM is empty. Run participant_build_geom first, e.g.\n' ...
         '  global WRIST_GEOM\n' ...
         '  WRIST_GEOM = participant_build_geom();']));
end
g = WRIST_GEOM;

sig = zeros(size(x));
epr = zeros(size(x));
[s0, e0] = tissueProps('soft', f);
sig(:) = s0;
epr(:) = e0;

env = g.envelope;

% ---- skin / fat / soft tissue rings -------------------------------
inSkin = inEllipse(x, y, env.cx, env.cy, env.aSkin,  env.bSkin,  env.phi) & ...
        ~inEllipse(x, y, env.cx, env.cy, env.aSkinIn, env.bSkinIn, env.phi);
inFat  = inEllipse(x, y, env.cx, env.cy, env.aSkinIn, env.bSkinIn, env.phi) & ...
        ~inEllipse(x, y, env.cx, env.cy, env.aFatIn,  env.bFatIn,  env.phi);

[s1, e1] = tissueProps('skin', f);
[s2, e2] = tissueProps('fat',  f);
sig(inSkin) = s1;  epr(inSkin) = e1;
sig(inFat)  = s2;  epr(inFat)  = e2;

% ---- bones ---------------------------------------------------------
[sb, eb] = tissueProps('bone', f);
for k = 1:numel(g.bonePolys)
    poly = g.bonePolys{k};
    in = inpolygon(x, y, poly(:, 1), poly(:, 2));
    sig(in) = sb;  epr(in) = eb;
end

% ---- cartilage (after bone, so bone wins on overlap) ---------------
[sc, ec] = tissueProps('cartilage', f);
for k = 1:numel(g.cartPolys)
    poly = g.cartPolys{k};
    in = inpolygon(x, y, poly(:, 1), poly(:, 2));
    sig(in) = sc;  epr(in) = ec;
end

% ---- median nerve --------------------------------------------------
[sn, en] = tissueProps('nerve', f);
if ~isempty(g.nerve)
    in = inEllipse(x, y, g.nerve(1), g.nerve(2), g.nerve(3), g.nerve(4), g.nerve(5));
    sig(in) = sn;  epr(in) = en;
end

% ---- flexor tendons -------------------------------------------------
[st, et] = tissueProps('tendon', f);
for k = 1:size(g.tendons, 1)
    in = ((x - g.tendons(k,1)).^2 + (y - g.tendons(k,2)).^2) <= g.tendons(k,3)^2;
    sig(in) = st;  epr(in) = et;
end

% ---- flexor retinaculum ----------------------------------------------
[sl, el] = tissueProps('ligament', f);
if size(g.retin, 1) >= 3
    in = inpolygon(x, y, g.retin(:, 1), g.retin(:, 2));
    sig(in) = sl;  epr(in) = el;
end
end


% ==================================================================== local ==
function tf = inEllipse(x, y, cx, cy, A, B, phi)
xr =  (x - cx)*cos(phi) + (y - cy)*sin(phi);
yr = -(x - cx)*sin(phi) + (y - cy)*cos(phi);
tf = (xr/A).^2 + (yr/B).^2 <= 1;
end


function [s, e] = tissueProps(name, f)
persistent T
if isempty(T)
    T.f = log10([1e2 1e3 1e4 1e5 1e6 1e7]);
    % rows: bone, soft, fat, skin, tendon, ligament, nerve, cartilage
    T.sigma = [ ...
        0.012 0.014 0.017 0.021 0.030 0.050;   % bone
        0.250 0.280 0.320 0.360 0.450 0.620;   % deep soft tissue
        0.020 0.021 0.023 0.025 0.030 0.045;   % fat
        0.020 0.050 0.120 0.250 0.420 0.550;   % skin
        0.320 0.360 0.400 0.440 0.520 0.600;   % tendon
        0.300 0.340 0.370 0.400 0.500 0.580;   % ligament
        0.170 0.220 0.260 0.300 0.420 0.550;   % nerve
        0.050 0.080 0.120 0.180 0.300 0.450];  % cartilage
    T.epsr = [ ...
        1e4   3e3   6e2   1.5e2 5e1   2e1;   % bone
        2e5   4e4   8e3   2.5e3 8e2   2.5e2; % deep soft tissue
        1e4   3e3   5e2   1e2   3e1   1.5e1; % fat
        1e5   2e4   4e3   1.2e3 4e2   1.5e2; % skin
        1.8e5 3.5e4 7e3   2.2e3 7e2   2.2e2; % tendon
        1.5e5 3e4   6e3   2e3   6e2   2e2;   % ligament
        2e5   4e4   8e3   2.5e3 8e2   2.5e2; % nerve
        1e4   4e3   1e3   3e2   1e2   4e1];  % cartilage
    T.names = {'bone','soft','fat','skin','tendon','ligament','nerve','cartilage'};
end
k = find(strcmp(T.names, name));
if isempty(k)
    error('wrist_tissue_sigma_real:tissue', 'Unknown tissue: %s', name);
end
lf = log10(f);
s = 10^interp1(T.f, log10(T.sigma(k, :)), lf, 'linear', 'extrap');
e = 10^interp1(T.f, log10(T.epsr(k, :)),  lf, 'linear', 'extrap');
end
