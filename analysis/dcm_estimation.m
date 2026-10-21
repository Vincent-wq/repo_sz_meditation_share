% =========================================================================
% SPECTRAL DYNAMIC CAUSAL MODELING (spDCM) BATCH SCRIPT
% =========================================================================
% This script loads the Python-exported .mat timecourses, constructs the 
% canonical SPM12 DCM struct for resting-state fMRI, and estimates the 
% effective (directed) connectivity matrices using Variational Laplace.

% 1. Setup Paths (UPDATE THESE TO YOUR LOCAL PATHS)
derive_dir = 'C:\PRIVATE_INPUTS\sz_meditation\derivatives\results_20260701';
data_dir = 'C:\PRIVATE_INPUTS\sz_meditation\derivatives\results_20260701\spDCM_Export_TR1p5_TE30'; 
out_dir  = 'C:\PRIVATE_INPUTS\sz_meditation\derivatives\results_20260701\spDCM_Estimated_TR1p5_TE30';
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

% Ensure SPM12 is initialized
spm('defaults', 'FMRI');
spm_jobman('initcfg');

% 2. Find all exported timecourse files
files = dir(fullfile(data_dir, '*_timecourses.mat'));
fprintf('Found %d files for spDCM estimation.\n', length(files));

% 3. Loop through each participant/session
for i = 1:length(files)
    
    file_name = files(i).name;
    file_path = fullfile(data_dir, file_name);
    fprintf('\nEstimating Model %d of %d: %s\n', i, length(files), file_name);
    
    % Load the Python data
    py_data = load(file_path);
    assert(isfield(py_data, 'TR') && isfield(py_data, 'TE'), ...
        'Re-export this input with both TR and TE.');
    TR = double(py_data.TR);
    TE = double(py_data.TE);
    validateattributes(TR, {'numeric'}, {'real','scalar','finite','positive'});
    validateattributes(TE, {'numeric'}, {'real','scalar','finite','positive'});

    Y_data  = py_data.Y;         % Timecourses (TRs x Networks)

    % MATLAB loads Python strings in a slightly complex cell format depending 
    % on the scipy.io version. This ensures names extract cleanly as a cell array of strings.
    if iscell(py_data.names)
        node_names = cellstr(char(py_data.names));
    else
        node_names = cellstr(py_data.names);
    end
    
    [v, n] = size(Y_data); % v = number of TRs, n = number of networks (9)

    % ==========================================
    % 4. BUILD THE DCM STRUCT
    % ==========================================
    clear DCM;
    
    % Data specifications
    DCM.Y.y    = Y_data;
    DCM.Y.dt   = TR;
    DCM.TE = TE;
    DCM.Y.name = node_names;
    DCM.v      = v;
    DCM.n      = n;

    DCM.acquisition.TR_seconds = TR;
    DCM.acquisition.TE_seconds = TE;
    out_name = fullfile(out_dir, ['DCM_', file_name]);
    assert(~exist(out_name, 'file'), 'Choose a fresh output path.');
    DCM.name = out_name;

    % Inputs (Empty for resting state)
    DCM.U.u    = zeros(v, 1);
    DCM.U.name = {'null'};
    
    % Architecture (The 'A' Matrix)
    % Setting ones(n,n) creates a "Fully Connected" model. 
    % We let SPM estimate all possible directed edges, and we will prune 
    % them later using Parametric Empirical Bayes (PEB).
    DCM.a = ones(n, n);      
    DCM.b = zeros(n, n, 0);  % No task modulations
    DCM.c = zeros(n, 1);     % No driving inputs
    DCM.d = zeros(n, n, 0);  % No non-linear modulations
    
    % spDCM Options
    DCM.options.analysis   = 'CSD';  % Trigger for Spectral DCM (Cross-Spectral Density)
    DCM.options.nonlinear  = 0;      
    DCM.options.two_state  = 0;
    DCM.options.stochastic = 0;
    DCM.options.centre     = 1;      % Mean-center time series
    DCM.options.maxnodes   = n;
    
    % Hemodynamic Priors
    DCM.delays = repmat(TR/2, n, 1); % Slice timing approximation
    
    % ==========================================
    % 5. ESTIMATE MODEL (Variational Laplace)
    % ==========================================
    % spm_dcm_fmri_csd computes the cross-spectra and fits the model
    try
        DCM = sz_spm_dcm_fmri_csd_TE(DCM);
        save(out_name, 'DCM');
        fprintf('Successfully estimated and saved %s\n', out_name);
    catch ME
        fprintf('FAILED to estimate %s\n', file_name);
        disp(ME.message);
    end
    
end

fprintf('\n=== ALL spDCM ESTIMATIONS COMPLETE ===\n');

%%
%save(fullfile(derive_dir, 'spDCM_Run_Workspace_TR1p5_TE30.mat'), '-v7.3');
