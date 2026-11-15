function rings = participant_cross_section(m, z0)
%PARTICIPANT_CROSS_SECTION  Intersect a surface mesh with plane z = z0 (mm).
%
%   rings = participant_cross_section(m, z0)
%
%   m     : struct from participant_read_inp (nodes 3xN, nid 1xN, elems PxK,
%           etype char). elems column 1 is the element ID; node IDs start at
%           column 2.
%   z0    : transverse plane position in mm (positive = distal / fingers)
%
%   rings : cell array of closed 2-D polygons (N x 2, x/y in mm)
%
%   Element types handled:
%     R3D3 / S3      : triangle shell            -> face [1 2 3]
%     C3D6 / C3D6H   : 6-node wedge (solid)      -> triangle cap [1 2 3]
%     M3D4R / S4     : quad shell                -> face [1 2 3 4]
%
%   For C3D6 cartilage wedges only the triangle cap is intersected on
%   purpose: the cartilage is a single wedge layer whose cap triangles form
%   one complete triangulated surface. Cutting the lateral quad faces as
%   well adds interior edges and splits the contour into many small rings.

switch upper(m.etype)
    case {'R3D3', 'S3', 'S3R', 'STRI65'}
        faces = {[1 2 3]};
    case {'C3D6', 'C3D6H', 'C3D6R'}
        faces = {[1 2 3]};
    case {'M3D4R', 'S4', 'S4R', 'S4R5'}
        faces = {[1 2 3 4]};
    otherwise
        faces = {[1 2 3]};
end

epsZ = 1e-8;
z0   = z0 + 1e-6;                 % tiny offset avoids vertices exactly on plane
segs = zeros(0, 4);
E = m.elems;

for k = 1:size(E, 1)
    eid = E(k, 2:end);             % column 1 is the element ID, rest are nodes
    eid = eid(eid > 0);
    for fi = 1:numel(faces)
        idx = faces{fi};
        idx = idx(idx <= numel(eid));
        if numel(idx) < 3
            continue;
        end
        [~, col] = ismember(eid(idx), m.nid);
        if any(col == 0)
            continue;
        end
        v = m.nodes(:, col);                 % 3 x nv
        zz = v(3, :);
        if all(zz < z0 - epsZ) || all(zz > z0 + epsZ)
            continue;
        end
        pts = zeros(2, 0);
        nv = numel(idx);
        for e = 1:nv
            a = e;
            b = mod(e, nv) + 1;
            za = zz(a);
            zb = zz(b);
            if (za < z0-epsZ && zb > z0+epsZ) || (za > z0+epsZ && zb < z0-epsZ)
                t = (z0 - za) / (zb - za);
                pts(:, end+1) = v(1:2, a)*(1-t) + v(1:2, b)*t;   %#ok<AGROW>
            elseif abs(za - z0) <= epsZ && (zb < z0-epsZ || zb > z0+epsZ)
                pts(:, end+1) = v(1:2, a);                       %#ok<AGROW>
            elseif abs(zb - z0) <= epsZ && (za < z0-epsZ || za > z0+epsZ)
                pts(:, end+1) = v(1:2, b);                       %#ok<AGROW>
            end
        end
        if size(pts, 2) == 2
            segs(end+1, :) = [pts(:,1).' pts(:,2).'];            %#ok<AGROW>
        elseif size(pts, 2) >= 3
            segs(end+1, :) = [pts(:,1).' pts(:,2).'];            %#ok<AGROW>
            segs(end+1, :) = [pts(:,1).' pts(:,3).'];            %#ok<AGROW>
        end
    end
end

rings = {};
if isempty(segs)
    return;
end

% ---------------- chain segments into closed loops ----------------
keyfun = @(p) sprintf('%.4f,%.4f', p(1), p(2));
nseg = size(segs, 1);
KA = cell(1, nseg);
KB = cell(1, nseg);
for i = 1:nseg
    KA{i} = keyfun(segs(i, 1:2));
    KB{i} = keyfun(segs(i, 3:4));
end

% deduplicate identical segments, drop zero-length ones
seen  = containers.Map('KeyType', 'char', 'ValueType', 'int32');
keep  = true(1, nseg);
for i = 1:nseg
    if strcmp(KA{i}, KB{i})
        keep(i) = false;
        continue;
    end
    s1 = sort({KA{i}, KB{i}});
    key = [s1{1} '|' s1{2}];
    if seen.isKey(key)
        keep(i) = false;
    else
        seen(key) = i;
    end
end
idxKeep = find(keep);
segs = segs(idxKeep, :);
KA   = KA(idxKeep);
KB   = KB(idxKeep);
nseg = numel(KA);

% adjacency: key -> unique neighbor keys
adj = containers.Map('KeyType', 'char', 'ValueType', 'any');
for i = 1:nseg
    if adj.isKey(KA{i})
        adj(KA{i}) = [adj(KA{i}), {KB{i}}];
    else
        adj(KA{i}) = {KB{i}};
    end
    if adj.isKey(KB{i})
        adj(KB{i}) = [adj(KB{i}), {KA{i}}];
    else
        adj(KB{i}) = {KA{i}};
    end
end
keys = adj.keys;
for i = 1:numel(keys)
    adj(keys{i}) = unique(adj(keys{i}));
end

used = false(1, nseg);
for i = 1:nseg
    if used(i)
        continue;
    end
    startK = KA{i};
    ringK  = {KA{i}, KB{i}};
    used(i) = true;
    cur  = KB{i};
    prev = KA{i};
    guard = 0;
    while ~strcmp(cur, startK) && guard < 2e6
        guard = guard + 1;
        nbrs = adj(cur);
        nxt = '';
        for j = 1:numel(nbrs)
            if ~strcmp(nbrs{j}, prev)
                nxt = nbrs{j};
                break;
            end
        end
        if isempty(nxt)
            break;
        end
        found = false;
        for j = 1:nseg
            if used(j)
                continue;
            end
            if (strcmp(KA{j}, cur) && strcmp(KB{j}, nxt)) || ...
               (strcmp(KB{j}, cur) && strcmp(KA{j}, nxt))
                used(j) = true;
                found = true;
                break;
            end
        end
        if ~found
            break;
        end
        ringK{end+1} = nxt;                     %#ok<AGROW>
        prev = cur;
        cur  = nxt;
    end
    if numel(ringK) >= 4
        XY = zeros(numel(ringK), 2);
        for j = 1:nseg
            mKa = strcmp(KA{j}, ringK);
            mKb = strcmp(KB{j}, ringK);
            if any(mKa)
                XY(mKa, :) = repmat(segs(j, 1:2), sum(mKa), 1);
            end
            if any(mKb)
                XY(mKb, :) = repmat(segs(j, 3:4), sum(mKb), 1);
            end
        end
        rings{end+1} = XY;                      %#ok<AGROW>
    end
end
end
