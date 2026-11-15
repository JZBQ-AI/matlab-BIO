function m = participant_read_inp(file)
%PARTICIPANT_READ_INP  Read Abaqus INP node/element blocks (R3D3/C3D6/M3D4R).
%
%   m = participant_read_inp(file)
%
%   m.nodes  : 3 x N double, (x, y, z) in mm
%   m.nid    : 1 x N double, node IDs
%   m.elems  : P x K double, element node IDs (padded with 0), first column is
%              the element ID
%   m.etype  : char, element type of the first *ELEMENT block
%   m.name   : char, base file name
%
%   Works for the PARTICIPANT_1_R dataset (R3D3 bone shells, C3D6 cartilage
%   wedges, M3D4R ligament membranes).

fid = fopen(file, 'r');
if fid < 0
    error('participant_read_inp:open', 'Cannot open file: %s', file);
end

nodes = zeros(3, 0);
nid   = zeros(1, 0);
elems = zeros(0, 1);
etype = '';
section = '';

while ~feof(fid)
    line = strtrim(fgetl(fid));
    if isempty(line) || strncmp(line, '**', 2)
        continue;
    end
    if line(1) == '*'
        section = '';
        if strncmp(line, '*NODE', 5)
            section = 'NODE';
        elseif contains(line, '*ELEMENT')
            section = 'ELEMENT';
            tok = regexp(line, 'TYPE=([A-Za-z0-9]+)', 'tokens', 'once');
            if ~isempty(tok)
                etype = tok{1};
            end
        end
        continue;
    end
    parts = strsplit(line, ',');
    if strcmp(section, 'NODE')
        v = str2double(parts);
        if numel(v) >= 4 && all(isfinite(v))
            nid(end+1) = v(1);                              %#ok<AGROW>
            nodes(:, end+1) = v(2:4).';                     %#ok<AGROW>
        end
    elseif strcmp(section, 'ELEMENT')
        v = str2double(parts);
        if numel(v) >= 2 && all(isfinite(v))
            row = zeros(1, numel(v));
            row(1:numel(v)) = v;
            if size(elems, 2) < numel(v)
                elems(:, end+1:numel(v)) = 0;
            end
            elems(end+1, 1:numel(v)) = row;                 %#ok<AGROW>
        end
    end
end
fclose(fid);

m.nodes = nodes;
m.nid   = nid;
m.elems = elems;
m.etype = etype;
[~, nm] = fileparts(file);
m.name  = nm;
end
