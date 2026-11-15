function geom = participant_build_geom(dataDir, z0, marginMm, includeSoft)
%PARTICIPANT_BUILD_GEOM  Build a 2-D wrist cross-section model from the
%   PARTICIPANT_1_R Abaqus dataset (bones, cartilage).
%
%   geom = participant_build_geom(dataDir, z0, marginMm, includeSoft)
%
%   dataDir     : folder that contains the PARTICIPANT_1_R sub-folders
%                 (BONE_GEOMETRY, CARTILAGE_GEOMETRY, ...)
%   z0          : transverse level in mm (default 10). Positive z = distal /
%                 toward the fingers. z = -10 gives the distal forearm
%                 (radius + ulna), z = 10..15 gives the carpal row.
%   marginMm    : soft-tissue + skin thickness added around the bone hull
%                 (default 11 mm)
%   includeSoft : true to add approximate flexor tendons, median nerve and
%                 flexor retinaculum (default true). PARTICIPANT_1_R does not
%                 contain tendon/nerve geometry, so these are placed from
%                 anatomical landmarks and are adjustable.
%
%   Returns geom with fields:
%     geom.z0, geom.margin, geom.unit ('mm')
%     geom.boneNames, geom.cartNames
%     geom.bonePolys : cell of Nx2 polygons (all bone rings, area filtered)
%     geom.cartPolys : same for cartilage
%     geom.envelope  : struct cx, cy, phi, aSkin, bSkin, aSkinIn, bSkinIn,
%                      aFatIn, bFatIn (half axes in mm)
%     geom.tendons   : 4 x 3 [cx cy r]
%     geom.nerve     : 1 x 5 [cx cy a b phi]
%     geom.retin     : N x 2 polygon
%     geom.gd, geom.sf : PDE geometry description matrix / set formula
%     geom.ns        : cell row vector of shape names matching gd columns
%                      (R2025b decsg requires a cellstr/string vector here)
%     geom.partNames : cell of names matching the gd columns
%     geom.repairLog : cellstr, rings that had to be rebuilt because they
%                      were self-intersecting after the plane cut
%
%   Geometry-description notes (why gd looks the way it does)
%   -------------------------------------------------------
%   decsg validates every column of gd (toolbox/shared/pdelib/decsg.m,
%   "validateGD"): the whole matrix must be finite, so NaN padding is rejected
%   outright by its mustBeFinite argument check; for a type-2 (polygon) column
%   rows 3..2+2*n hold the vertices and every row after that must be exactly
%   zero; and each polygon must be a simple closed curve with unique vertices
%   (Shamos-Hoey sweep). The mesh cross-sections do violate these rules at
%   grazing cuts, so every ring is passed through sanitizePoly and the gd
%   matrix is zero-padded instead of being resampled to a fixed vertex count.

if nargin < 1 || isempty(dataDir)
    dataDir = defaultDataDir();
end
if ~exist(fullfile(dataDir, 'BONE_GEOMETRY'), 'dir')
    error('participant_build_geom:dir', ...
        sprintf(['PARTICIPANT_1_R not found automatically. Pass the folder explicitly:\n' ...
         '  geom = participant_build_geom(''<path>\\PARTICIPANT_1_R'');']));
end
if nargin < 2 || isempty(z0),        z0        = 10; end
if nargin < 3 || isempty(marginMm),  marginMm  = 11; end
if nargin < 4 || isempty(includeSoft), includeSoft = true; end

geom.z0      = z0;
geom.margin  = marginMm;
geom.unit    = 'mm';
geom.dataDir = dataDir;
geom.repairLog = {};

% ------------------------------------------------------------------ bones --
boneDir = fullfile(dataDir, 'BONE_GEOMETRY');
cartDir = fullfile(dataDir, 'CARTILAGE_GEOMETRY');
if ~exist(boneDir, 'dir') || ~exist(cartDir, 'dir')
    error('participant_build_geom:dir', ...
        sprintf(['BONE_GEOMETRY/CARTILAGE_GEOMETRY not found under:\n  %s\n' ...
         'Pass the PARTICIPANT_1_R folder explicitly via ''dataDir''.'], dataDir));
end

[bonePolys, boneNames, logB] = loadCrossPolys(boneDir, z0, 'BONE');
[cartPolys, cartNames, logC] = loadCrossPolys(cartDir, z0, 'CARTILAGE');

geom.boneNames = boneNames;
geom.cartNames = cartNames;
geom.bonePolys = bonePolys;
geom.cartPolys = cartPolys;
geom.repairLog = [logB; logC];

if isempty(bonePolys)
    error('participant_build_geom:nosection', ...
        'No bone cross-section found at z0 = %.1f mm. Choose another level.', z0);
end

% ----------------------------------------------------------------- envelope --
geom.envelope = fitEnvelope(bonePolys, cartPolys, marginMm);

% ---------------------------------------------------- soft structures -------
geom.tendons = zeros(0, 3);
geom.nerve   = zeros(1, 5);
geom.retin   = zeros(0, 2);
if includeSoft
    [geom.tendons, geom.nerve, geom.retin] = placeSoft(geom);
end

% ------------------------------------------------------------ PDE geometry --
[gd, ns, sf, partNames] = buildPDE(geom);
geom.gd        = gd;
geom.ns        = ns;
geom.sf        = sf;
geom.partNames = partNames;

% quick stats
allx = []; ally = [];
for k = 1:numel(bonePolys)
    allx = [allx; bonePolys{k}(:,1)]; %#ok<AGROW>
    ally = [ally; bonePolys{k}(:,2)]; %#ok<AGROW>
end
geom.stats.boneSpanX_mm = max(allx) - min(allx);
geom.stats.boneSpanY_mm = max(ally) - min(ally);
geom.stats.skinA_mm     = geom.envelope.aSkin;
geom.stats.skinB_mm     = geom.envelope.bSkin;
geom.stats.nBoneRegions = numel(bonePolys);
geom.stats.nCartRegions = numel(cartPolys);
end


% ==================================================================== local ==
function [polys, names, log] = loadCrossPolys(folder, z0, tag)
files = dir(fullfile(folder, '*.INP'));
polys = cell(0, 1);
names = cell(0, 1);
log   = cell(0, 1);
for fi = 1:numel(files)
    m = participant_read_inp(fullfile(folder, files(fi).name));
    if isempty(m.nodes) || isempty(m.elems)
        continue;
    end
    zz = m.nodes(3, :);
    if max(zz) < z0 - 1 || min(zz) > z0 + 1
        continue;                          % bone/cartilage outside this level
    end
    rings = participant_cross_section(m, z0);
    base = regexprep(files(fi).name, 'PARTICIPANT_1_R_|_GEOM\.INP|\.INP', '', 'ignorecase');
    for r = 1:numel(rings)
        nm = sprintf('%s_%s', tag, base);
        [poly, why] = sanitizePoly(rings{r});
        rawArea = 0;
        if size(rings{r}, 1) >= 3
            rawArea = abs(polyArea(rings{r}(:, 1:2)));
        end
        if size(poly, 1) >= 3 && abs(polyArea(poly)) >= 1.0
            polys{end+1, 1} = poly;   %#ok<AGROW>
            names{end+1, 1} = nm;     %#ok<AGROW>
            if ~isempty(why)
                log{end+1, 1} = sprintf('%s[%d] %s', nm, r, why); %#ok<AGROW>
            end
        elseif rawArea >= 1.0
            % a region that should have survived was lost: make it visible
            log{end+1, 1} = sprintf('%s[%d] DROPPED (raw area %.2f mm^2): %s', ...
                nm, r, rawArea, why); %#ok<AGROW>
        end
    end
end
end


% ------------------------------------------------------------- polygon prep --
function [poly, why] = sanitizePoly(ring)
%SANITIZEPOLY  Turn one raw plane-cut ring into a strictly simple polygon.
%   decsg requires every type-2 column to be a simple closed curve with
%   unique vertices. Raw mesh contours break that in three ways: repeated
%   vertices (shared cut points), zero-length/duplicate edges, and - at
%   grazing cuts - genuinely self-crossing loops produced when the segment
%   chaining has to choose between several branches at a junction.
%   Returns poly (possibly empty when the ring is degenerate) and a short
%   text reason when the ring had to be repaired or dropped.

poly = zeros(0, 2);
why  = '';
if isempty(ring), return; end
P = ring(:, 1:2);
P = P(all(isfinite(P), 2), :);
if size(P, 1) < 3, return; end

aScale = max(1, max(abs(P(:))));
% tolDup: raw contours contain sub-micron vertex spacings (0.3-0.8 um) left by
% the plane cut; the R2025b mesher rejects the degenerate edges they create
% ("Geometry error"), so vertices closer than ~10-30 um are merged. That is
% 1-3% of the default Hmax (0.8 mm) and cannot change the anatomy.
tolDup = max(1e-2, 1e-3 * aScale);
% tolTst: only near-exact coincidences matter for the self-intersection test.
tolTst = 1e-6 * aScale;

P = dropConsecutiveDuplicates(P, tolDup);
P = dropCollinear(P, tolDup);
if size(P, 1) < 3, return; end

[nCross, nTouch] = polyViolations(P, tolTst);
if nCross == 0 && nTouch == 0
    poly = orientCCW(P);
    return;
end

% The ring crosses (or touches) itself: rebuild it from the largest simple
% region of its polyshape decomposition. polyshape resolves the crossings and
% returns regions whose boundaries are simple by construction, so its output
% is trusted without re-running the (much stricter) local test above.
a0 = abs(polyArea(P));
Q  = rebuildSimple(P);
if isempty(Q)
    why = sprintf('dropped: %d crossings / %d touches, no simple region left (area %.3f)', ...
        nCross, nTouch, a0);
    return;
end
Q = dropConsecutiveDuplicates(Q, tolDup);
Q = dropCollinear(Q, tolDup);
if size(Q, 1) < 3
    why = sprintf('dropped: repaired ring degenerate (area %.3f)', a0);
    return;
end
a1 = abs(polyArea(Q));
if a1 < 0.5 * a0
    why = sprintf('repaired: bow-tie ring, area %.3f -> %.3f mm^2', a0, a1);
else
    why = sprintf('repaired: %d crossings removed, area %.3f -> %.3f mm^2', nCross, a0, a1);
end
poly = orientCCW(Q);
end


function P = dropConsecutiveDuplicates(P, tol)
if size(P, 1) < 2, return; end
guard = 0;
while guard < 10 * size(P, 1) + 10
    guard = guard + 1;
    if size(P, 1) < 4, break; end
    d = sqrt(sum(diff(P, 1, 1).^2, 2));
    if all(d > tol)
        break;
    end
    keep = [true; d > tol];
    P = P(keep, :);
end
% the closing vertex must not repeat the first one
guard = 0;
while size(P, 1) > 3 && norm(P(1, :) - P(end, :)) <= tol && guard < 100
    P(end, :) = [];
    guard = guard + 1;
end
end


function P = dropCollinear(P, tol)
%DROPCOLLINEAR  Remove vertices that lie on the straight line between their
%   neighbours. The outline is unchanged, it only removes zero-area corners
%   (and the degenerate edges they create) before decsg sees the polygon.
if size(P, 1) < 4, return; end
i = 1;
guard = 0;
maxIter = 10 * size(P, 1) + 10;
while guard < maxIter
    guard = guard + 1;
    n = size(P, 1);
    if n < 4, break; end
    i  = mod(i - 1, n) + 1;
    ip = mod(i - 2, n) + 1;
    in = mod(i, n) + 1;
    a = P(ip, :);  b = P(i, :);  c = P(in, :);
    v = c - a;
    L = norm(v);
    if L <= tol                       % neighbours coincide: b is redundant
        P(i, :) = [];
        continue;
    end
    d = abs((b(1) - a(1)) * v(2) - (b(2) - a(2)) * v(1)) / L;    % distance to line a-c
    t = ((b(1) - a(1)) * v(1) + (b(2) - a(2)) * v(2)) / (L * L); % projection onto a-c
    if d <= tol && t > 0 && t < 1
        P(i, :) = [];                 % b lies on the segment: drop it
    else
        i = i + 1;
    end
end
end


function [nCross, nTouch] = polyViolations(P, tol)
%POLYVIOLATIONS  Count proper self-crossings and self-touching edge pairs.
n = size(P, 1);
nCross = 0;
nTouch = 0;
if n < 4, return; end
for i = 1:n
    a1 = P(i, :);
    a2 = P(mod(i, n) + 1, :);
    for j = i + 1:n
        if mod(j, n) + 1 == i || j == mod(i, n) + 1
            continue;                      % adjacent edges share a vertex
        end
        b1 = P(j, :);
        b2 = P(mod(j, n) + 1, :);
        d1 = cross2(a1, a2, b1);
        d2 = cross2(a1, a2, b2);
        d3 = cross2(b1, b2, a1);
        d4 = cross2(b1, b2, a2);
        if ((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) && ...
           ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0))
            nCross = nCross + 1;
        elseif (abs(d1) <= tol && onSeg(a1, a2, b1, tol)) || ...
               (abs(d2) <= tol && onSeg(a1, a2, b2, tol)) || ...
               (abs(d3) <= tol && onSeg(b1, b2, a1, tol)) || ...
               (abs(d4) <= tol && onSeg(b1, b2, a2, tol))
            nTouch = nTouch + 1;
        end
    end
end
end


function d = cross2(p, q, r)
d = (q(1) - p(1)) * (r(2) - p(2)) - (q(2) - p(2)) * (r(1) - p(1));
end


function tf = onSeg(a, b, p, tol)
%ONS  Is p (already known to be on the line through a,b) between a and b?
d = b - a;
L2 = d(1)^2 + d(2)^2;
if L2 <= 0, tf = false; return; end
t = ((p(1) - a(1)) * d(1) + (p(2) - a(2)) * d(2)) / L2;
tolT = tol / sqrt(L2);
tf = t > tolT && t < 1 - tolT;
end


function Q = rebuildSimple(P)
%REBUILDSIMPLE  Largest simple region of a self-intersecting ring.
Q = zeros(0, 2);
if exist('polyshape', 'file') ~= 2
    return;                                % no polyshape: drop the ring
end
% polyshape warns when its input self-intersects; that is exactly the case
% this function exists for, so keep the command window quiet.
ws = warning('off', 'all');
try
    pg = polyshape(P(:, 1), P(:, 2), 'Simplify', true, 'KeepCollinearPoints', false);
    if pg.NumRegions == 0
        warning(ws);
        return;
    end
    R = regions(pg);
    if isempty(R)
        warning(ws);
        return;
    end
    [aMax, k] = max(area(R));
    if ~isfinite(aMax) || aMax <= 0
        warning(ws);
        return;
    end
    Q = R(k).Vertices;
catch
    Q = zeros(0, 2);
end
warning(ws);
if size(Q, 1) < 4
    Q = zeros(0, 2);
    return;
end
guard = 0;
while size(Q, 1) > 3 && norm(Q(1, :) - Q(end, :)) <= 1e-9 && guard < 10
    Q(end, :) = [];
    guard = guard + 1;
end
end


function P = orientCCW(P)
if polyArea(P) < 0
    P = flipud(P);
end
end


function a = polyArea(poly)
x = poly(:, 1);
y = poly(:, 2);
a = 0.5 * sum(x .* circshift(y, -1) - circshift(x, -1) .* y);
end


function env = fitEnvelope(bonePolys, cartPolys, marginMm)
P = [];
for k = 1:numel(bonePolys)
    P = [P; bonePolys{k}]; %#ok<AGROW>
end
for k = 1:numel(cartPolys)
    P = [P; cartPolys{k}]; %#ok<AGROW>
end
mu = mean(P, 1);
Pc = P - mu;
[V, D] = eig(cov(Pc));
[~, order] = sort(diag(D), 'descend');
V = V(:, order);                       % principal axes
Q = Pc * V;                            % projections
amin = min(Q, [], 1);
amax = max(Q, [], 1);
a = 0.5 * (amax - amin);               % half axes along principal dirs
centerQ = 0.5 * (amax + amin);
center = (mu + (V * centerQ.').');           % world center of the fitted ellipse
phi = atan2(V(2, 1), V(1, 1));         % rotation of first axis

env.cx     = center(1);
env.cy     = center(2);
env.phi    = phi;
env.aSkin  = a(1) + marginMm;
env.bSkin  = a(2) + marginMm;
env.aSkinIn = 0.94 * env.aSkin;
env.bSkinIn = 0.94 * env.bSkin;
env.aFatIn  = 0.87 * env.aSkin;
env.bFatIn  = 0.87 * env.bSkin;
end


function [tendons, nerve, retin] = placeSoft(geom)
P = [];
for k = 1:numel(geom.bonePolys)
    P = [P; geom.bonePolys{k}]; %#ok<AGROW>
end
for k = 1:numel(geom.cartPolys)
    P = [P; geom.cartPolys{k}]; %#ok<AGROW>
end
xmin = min(P(:,1)); xmax = max(P(:,1));
ymin = min(P(:,2)); ymax = max(P(:,2));
w = xmax - xmin;
h = ymax - ymin;

% assume +y = palmar (volar) side; place structures in the palmar corridor
xLo = xmin + 0.30*w;
xHi = xmax - 0.25*w;
yT  = ymin + 0.72*h;

xs = linspace(xLo, xHi, 4);
r  = min(2.8, 0.09*w);
tendons = [xs(:), repmat(yT, 4, 1), repmat(r, 4, 1)];

nerve = [xs(1) - 3.5, yT - 1.5, 3.2, 1.8, 0];

% thin band of the flexor retinaculum across the palmar opening
retin = [xLo - 2, yT + 4;
         xHi + 2, yT + 4;
         xHi + 2, yT + 6;
         xLo - 2, yT + 6];
end


function [gd, ns, sf, partNames] = buildPDE(geom)
%BUILDPDE  Geometry description matrix for decsg / PDE Modeler.
%
% decsg (toolbox/shared/pdelib/decsg.m) rules that used to be broken here:
%   * mustBeFinite on gd  -> NaN padding (the classic mixed-length trick) is
%     rejected, so the matrix is zero-padded instead. decsg only reads rows
%     3..2+2*n of a type-2 column and requires the remaining rows to be zero.
%   * Shamos-Hoey self-intersection test per column -> every polygon must be
%     simple. Resampling all shapes to 128 arc-length points (the previous
%     workaround) *propagated* self-crossings from the raw cut instead of
%     removing them, and sampled chords could cross thin concave necks, so
%     the polygons now keep their true vertices after sanitizePoly.
Nv = 128;                                    % vertices for the ellipse approximations
cols = {};
partNames = {};

% envelope ellipses, approximated by Nv-gons (keeps the outer boundary a
% polygon so the electrode detection thresholds stay meaningful)
[cols, partNames] = addEllipse(cols, partNames, ...
    geom.envelope.cx, geom.envelope.cy, geom.envelope.aSkin, geom.envelope.bSkin, ...
    geom.envelope.phi, Nv, 'SKO');
[cols, partNames] = addEllipse(cols, partNames, ...
    geom.envelope.cx, geom.envelope.cy, geom.envelope.aSkinIn, geom.envelope.bSkinIn, ...
    geom.envelope.phi, Nv, 'SKI');
[cols, partNames] = addEllipse(cols, partNames, ...
    geom.envelope.cx, geom.envelope.cy, geom.envelope.aFatIn, geom.envelope.bFatIn, ...
    geom.envelope.phi, Nv, 'FAT');

% bones: true cross-section vertices
for k = 1:numel(geom.bonePolys)
    [cols, partNames] = addPoly(cols, partNames, geom.bonePolys{k}, sprintf('B%02d', k));
end

% cartilage: true cross-section vertices
for k = 1:numel(geom.cartPolys)
    [cols, partNames] = addPoly(cols, partNames, geom.cartPolys{k}, sprintf('C%02d', k));
end

% tendons (circles expressed as Nv-gons, a == b and phi == 0)
for k = 1:size(geom.tendons, 1)
    [cols, partNames] = addEllipse(cols, partNames, ...
        geom.tendons(k, 1), geom.tendons(k, 2), geom.tendons(k, 3), geom.tendons(k, 3), ...
        0, Nv, sprintf('T%02d', k));
end

% median nerve (ellipse expressed as polygon)
if numel(geom.nerve) >= 5
    [cols, partNames] = addEllipse(cols, partNames, ...
        geom.nerve(1), geom.nerve(2), geom.nerve(3), geom.nerve(4), ...
        geom.nerve(5), Nv, 'NERVE');
end

% flexor retinaculum (axis-aligned rectangle, kept as its 4 vertices)
if size(geom.retin, 1) >= 4
    [cols, partNames] = addPoly(cols, partNames, geom.retin, 'RETN');
end

if isempty(cols)
    error('participant_build_geom:empty', 'No valid geometry could be built.');
end

% zero-pad every column to the longest one (decsg reads only rows 3..2+2n)
nRows = max(cellfun(@numel, cols));
gd = zeros(nRows, numel(cols));
for k = 1:numel(cols)
    c = cols{k};
    gd(1:numel(c), k) = c;
end
if any(~isfinite(gd(:)))
    error('participant_build_geom:nonfinite', 'Geometry matrix contains non-finite values.');
end
% decsg needs row 2 == number of vertices and zeros beyond the last vertex
for k = 1:size(gd, 2)
    n = gd(2, k);
    if 2 + 2*n > nRows || any(gd(3 + 2*n:end, k))
        error('participant_build_geom:gdformat', ...
            'Malformed geometry column %d (%s).', k, partNames{k});
    end
end
sf = strjoin(partNames, '+');
% NS: R2025b decsg only accepts a string/cellstr vector here (a char matrix
% with one name per ROW is rejected by validateSFNS with GDNSMismatch), so the
% name space is returned as a cell row vector. decsg(gd, sf, ns) takes it as is.
ns = partNames;
end


function [cols, names] = addEllipse(cols, names, xc, yc, a, b, phi, N, name)
th = linspace(0, 2*pi, N+1);
th(end) = [];
xr = a*cos(th);
yr = b*sin(th);
x = xc + xr*cos(phi) - yr*sin(phi);
y = yc + xr*sin(phi) + yr*cos(phi);
xy = [x(:), y(:)];
cols{end+1} = [2; N; xy(:, 1); xy(:, 2)];
names{end+1} = name;
end


function [cols, names] = addPoly(cols, names, poly, name)
%ADDPOLY  Append one simple polygon as a type-2 geometry column.
%   Keeps the true vertices: decsg does not require a fixed vertex count, so
%   there is no decimation and no arc-length resampling (both of which used
%   to turn thin concave rings into self-intersecting ones).
if size(poly, 1) < 3
    return;
end
cols{end+1} = [2; size(poly, 1); poly(:, 1); poly(:, 2)];
names{end+1} = name;
end


function d = defaultDataDir()
% defaultDataDir  Auto-detect the PARTICIPANT_1_R folder (name is ASCII, so
% the Chinese parent folder can be found without hard-coding it).
d = '';
roots = {'E:\BaiduNetdiskDownload', 'E:\', 'D:\', 'F:\', ...
         'C:\Users\Administrator\Downloads'};
for k = 1:numel(roots)
    if exist(roots{k}, 'dir')
        d = findParticipant(roots{k});
        if ~isempty(d)
            return;
        end
    end
end
end


function d = findParticipant(root)
d = '';
items = dir(root);
for k = 1:numel(items)
    if items(k).isdir && strcmp(items(k).name, 'PARTICIPANT_1_R')
        d = fullfile(root, items(k).name);
        return;
    end
end
% one level deeper (e.g. E:\BaiduNetdiskDownload\<chinese>\PARTICIPANT_1_R)
for k = 1:numel(items)
    if items(k).isdir && ~strncmp(items(k).name, '.', 1)
        sub = fullfile(root, items(k).name);
        items2 = dir(sub);
        for j = 1:numel(items2)
            if items2(j).isdir && strcmp(items2(j).name, 'PARTICIPANT_1_R')
                d = fullfile(sub, items2(j).name);
                return;
            end
        end
    end
end
end
