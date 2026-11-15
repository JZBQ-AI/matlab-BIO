function run_wrist_real(varargin)
%RUN_WRIST_REAL  One-shot runner for the PARTICIPANT_1_R real-geometry model.
%
%   run_wrist_real                     % preview levels + sweep at z0 = 10
%   run_wrist_real('z0', -10)          % distal forearm level
%   run_wrist_real('freqs', logspace(2,5,13))
%
%   Step 1: build geometry (auto-detect the PARTICIPANT_1_R folder).
%   Step 2: preview cross-sections at several levels.
%   Step 3: frequency sweep with the real bone/cartilage geometry.
%   Step 4: PDE Modeler GUI script (jianmo_real) is available for manual work.

opt = struct('z0', 10, 'freqs', logspace(2, 6, 25));
for k = 1:2:numel(varargin)
    if ~isfield(opt, varargin{k})
        error('run_wrist_real:badOption', 'Unknown option: %s', num2str(varargin{k}));
    end
    opt.(varargin{k}) = varargin{k+1};
end

fprintf('\n===== PARTICIPANT_1_R real wrist model =====\n');
global WRIST_GEOM
WRIST_GEOM = participant_build_geom();      % auto-detect data folder

fprintf('\n--- preview cross-sections ---\n');
wrist_cross_preview();

fprintf('\n--- frequency sweep (z0 = %.0f mm) ---\n', opt.z0);
wrist_z_spectrum_real('z0', opt.z0, 'freqs', opt.freqs);

fprintf('\nDone. For the PDE Modeler GUI version run:  jianmo_real\n');
end
