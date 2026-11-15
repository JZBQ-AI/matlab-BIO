function c = wrist_sigma_coef_real(p, t, u, time)
%WRIST_SIGMA_COEF_REAL  Coefficient function for the PDE Modeler GUI.
%
%   Called by PDE Toolbox as c = f(p, t, u, time). Returns per-element
%   conductivity (real, S/m). For a complex solve (with permittivity) use
%   wrist_z_spectrum_real instead.
%
%   Frequency is read from the global WRIST_FREQ_HZ (default 100 kHz).
%   Geometry is read from the global WRIST_GEOM (participant_build_geom).

global WRIST_FREQ_HZ WRIST_GEOM
if isempty(WRIST_FREQ_HZ)
    WRIST_FREQ_HZ = 1e5;
end
if isempty(WRIST_GEOM)
    WRIST_GEOM = participant_build_geom();
end

nt = size(t, 2);
xc = mean(reshape(p(1, t(1:3, :)), 3, nt), 1);
yc = mean(reshape(p(2, t(1:3, :)), 3, nt), 1);
[sig, ~] = wrist_tissue_sigma_real(xc, yc, WRIST_FREQ_HZ);
c = sig;
end
