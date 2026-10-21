function peb_section3_final
% Fit the Section 3 longitudinal PEB models and their planned sensitivities.
% Input: estimated participant DCMs, group metadata and visit-specific mean FD.
% Output: one final model MAT and four scientific CSVs in data/dcm/PEB_section3.

%% 1. Input files and reporting settings
% Use the existing DCM estimates; TR and TE are read from them, not changed here.
project_dir = 'C:\PRIVATE_INPUTS\repo_sz_meditation';
dcm_dir = ['C:\PRIVATE_INPUTS\sz_meditation\derivatives\' ...
    'results_20260701\spDCM_Estimated_TR1p5_TE30'];
data_dir = fullfile(project_dir, 'data');
clinical_file = fullfile(data_dir, 'sz_meditation_all_DKT_20240724.csv');
motion_file = fullfile(data_dir, 'IC_Stats_2026.csv');
output_dir = fullfile(data_dir, 'dcm', 'PEB_section3');
%%
fit_models = true; % Set false to rebuild the four CSVs from the final model MAT.
model_file = fullfile(output_dir, 'PEB_models.mat');
if ~isfolder(output_dir), mkdir(output_dir); end

spm('defaults', 'FMRI');
spm_get_defaults('dcm.verbose', true);

if fit_models
    %% 2. Match participants, visits, treatment group and motion
    % Keep all available fitted scans; the main PEB uses all three fitted visits.
    clinical = readtable(clinical_file, 'TextType', 'string');
    motion = readtable(motion_file, 'TextType', 'string');
    subjects = unique(strtrim(string(clinical.participant_id)));
    n_subjects = numel(subjects);
    months = [0;3;8];
    clinical_ids = strtrim(string(clinical.participant_id));
    motion_ids = strtrim(string(motion.participant_id));
    group_map = containers.Map({'control','meditation'}, {-1,1});
    groups = nan(n_subjects,1);
    baseline_record = false(n_subjects,1);
    FD = nan(n_subjects,3);
    files = cell(n_subjects,3);
    DCMs = cell(n_subjects,3);
    for i = 1:n_subjects
        base = clinical_ids == subjects(i) & clinical.session == 0;
        baseline_record(i) = any(base);
        arm = unique(lower(strtrim(string(clinical.group(clinical_ids == subjects(i))))));
        arm = arm(~ismissing(arm) & arm ~= "");
        groups(i) = group_map(char(arm));
        for j = 1:3
            fd_row = motion_ids == subjects(i) & motion.session == months(j);
            if any(fd_row), FD(i,j) = motion.mean_fd(fd_row); end
            filename = fullfile(dcm_dir, sprintf('DCM_%s_ses-%d_timecourses.mat', ...
                char(subjects(i)), months(j)));
            if ~isfile(filename), continue; end
            files{i,j} = filename;
            loaded = load(filename, 'DCM');
            DCMs{i,j} = dcm_A_summary(loaded.DCM, filename);
        end
    end
    available = ~cellfun(@isempty, DCMs);
    discovery_index = find(all(available,2) & baseline_record);
    n_reference = numel(discovery_index);
    group_reference = groups(discovery_index);

    PEB = struct;
    PEB.Node_IDs = {'IC00_SMN'; 'IC02_SalVentA'; 'IC04_Vis_V1'; ...
        'IC05_aDMN'; 'IC09_pDMN'; 'IC15_DAN'; 'IC16_FPCN'; ...
        'IC18_SalVentB'; 'IC19_ECN'};
    PEB.Node_labels = {'SMN'; 'SMG–VAN'; 'mVIS'; 'mPFC–DMN'; 'PCu–DMN'; ...
        'SPL–DAN'; 'rIPL–FPN'; 'MFG–VAN/FPN'; 'lIFG–FPN'};
    PEB.Subjects = cellstr(subjects);
    PEB.Groups = groups;
    PEB.FD = FD;
    PEB.Available_scans = available;
    PEB.Discovery_index = discovery_index;
    PEB.DCM_files = files;
    PEB.DCM_summaries = DCMs;
    PEB.Time = months;
    PEB.Threshold = 0.95;
    PEB.BMA_seed = 24;
    PEB.Group_center = mean(group_reference);
    PEB.Inputs = struct('DCM_directory', dcm_dir, ...
        'Clinical_CSV', clinical_file, 'Motion_CSV', motion_file);
    PEB.SPM_version = spm('Ver');
    PEB.SPM_functions = struct('PEB', which('spm_dcm_peb'), ...
        'BMC', which('spm_dcm_peb_bmc'), 'BMR', which('spm_dcm_bmr_all'));
    PEB.Reduction = struct('NoFD', 'spm_dcm_peb_bmc(PEB)', ...
        'FD', 'spm_dcm_bmr_all: Mean/Group reduced; FD retained');

    %% 3. Fit the three-visit hierarchy without motion covariates
    % Center linear time and orthogonalize quadratic time; retain all A parameters.
    linear_time = months - mean(months);
    quadratic_time = linear_time.^2 - linear_time * ...
        ((linear_time'*linear_time) \ (linear_time'*(linear_time.^2)));
    quadratic_time = quadratic_time - mean(quadratic_time);
    Xwithin = [ones(3,1), linear_time, quadratic_time];
    Xbetween = [ones(n_reference,1), group_reference-PEB.Group_center];
    PEB.Xwithin = Xwithin;
    PEB.Xbetween_NoFD = Xbetween;

    % SPM uses a temporary file internally; keep it outside the data folder.
    original_dir = pwd;
    work_dir = tempname;
    mkdir(work_dir);
    restore_dir = onCleanup(@() finish_spm_work(original_dir, work_dir));
    cd(work_dir);
    PEB.Within = cell(n_reference,1);
    for i = 1:n_reference
        scans = files(discovery_index(i),:)';
        PEB.Within{i} = spm_dcm_peb(scans, Xwithin, {'A'});
    end
    M = struct('X', Xbetween, 'Xnames', {{'Mean','Group'}});
    PEB.NoFD = spm_dcm_peb(PEB.Within, M);
    rng(PEB.BMA_seed, 'twister');
    PEB.BMA_NoFD = spm_dcm_peb_bmc(PEB.NoFD);

    %% 4. Fit motion-adjusted PEB and retain motion terms during reduction
    % Reuse the same within-person PEBs and participants; FD is scaled by visit.
    motion_reference = FD(discovery_index,:);
    PEB.FD_center = mean(motion_reference,1);
    PEB.FD_scale = std(motion_reference,1,1);
    motion_z = (motion_reference-PEB.FD_center)./PEB.FD_scale;
    PEB.Xbetween_FD = [Xbetween, motion_z];
    M = struct('X', PEB.Xbetween_FD, ...
        'Xnames', {{'Mean','Group','FD0','FD3','FD8'}}, 'maxit', 256);
    PEB.FullFD = spm_dcm_peb(PEB.Within, M);
    packed = PEB.FullFD;
    packed.Ep = struct('Target', full(PEB.FullFD.Ep(:,1:2)), ...
        'FD', full(PEB.FullFD.Ep(:,3:5)));
    packed.M.pE = spm_unvec(spm_vec(PEB.FullFD.M.pE), packed.Ep);
    rng(PEB.BMA_seed, 'twister');
    averaged = spm_dcm_bmr_all(packed, {'Target'}, 'BMA');
    PEB.BMA_FD = averaged;
    PEB.BMA_FD.Ep = reshape(spm_vec(averaged.Ep), size(PEB.FullFD.Ep));
    PEB.BMA_FD.Pp = reshape(spm_vec(averaged.Pp), size(PEB.FullFD.Ep));
    PEB.BMA_FD.M.pE = PEB.FullFD.M.pE;

    %% 5. Fit window-specific motion and paired-visit sensitivities
    % Keep the mapped priors and discovery FD scale when including extra pairs.
    PEB.Paired = fit_paired_sensitivities(DCMs, subjects, groups, FD, ...
        discovery_index, Xwithin, PEB.Within, PEB.Node_IDs);
    clear restore_dir

    %% 6. Save the final fitted models for later inference
    % Preserve full posterior covariance and the original A marginals, not time series.
    PEB.Estimation_note = ['Native SPM estimates; code consolidation does not ' ...
        'certify numerical stationarity or remove the historical fit qualifications.'];
    save(model_file, 'PEB', '-v7');
else
    saved = load(model_file, 'PEB');
    PEB = saved.PEB;
end

%% 7. Export scientific parameters, contrasts and individual ECs
% Full-PEB intervals stay with full-PEB means; BMA means and selection stay separate.
[coefficients, contrasts, edge_dictionary, scan_table] = make_peb_tables(PEB);
writetable(coefficients, fullfile(output_dir, 'PEB_treatment_coefficients.csv'));
writetable(contrasts, fullfile(output_dir, 'PEB_contrasts.csv'));
writetable(edge_dictionary, fullfile(output_dir, 'Selected_edge_dictionary.csv'));
writetable(scan_table, fullfile(output_dir, 'DCM_selected_edges_all_available.csv'));
fprintf('\nPEB discovery: N=%d (%d control, %d meditation).\n', ...
    numel(PEB.Discovery_index), sum(PEB.Groups(PEB.Discovery_index)==-1), ...
    sum(PEB.Groups(PEB.Discovery_index)==1));
fprintf('Final models and four reporting tables: %s\n', output_dir);
end

%% Functions defined
function summary = dcm_A_summary(DCM, filename)
% Retain original A posterior/prior marginals for paired PEB and future inference.
    indices = spm_find_pC(DCM, {'A'});
    if isstruct(DCM.M.pC)
        prior_covariance = diag(spm_vec(DCM.M.pC));
    else
        prior_covariance = DCM.M.pC;
    end
    summary = struct('name', filename, 'Ep', struct('A', DCM.Ep.A), ...
        'Cp', full(DCM.Cp(indices,indices)), 'F', DCM.F);
    summary.M = struct('pE', struct('A', DCM.M.pE.A), ...
        'pC', full(prior_covariance(indices,indices)));
    summary.A_indices_in_original_DCM = indices;
    summary.Node_names = DCM.Y.name;
    summary.Acquisition = struct('TR_seconds', DCM.Y.dt, ...
        'TE_seconds', DCM.TE, 'Model_TE_seconds', DCM.M.TE);
end


function finish_spm_work(original_dir, work_dir)
% Restore the working directory and remove only this run's native temporary file.
    cd(original_dir);
    native_file = fullfile(work_dir, 'tmp.mat');
    if isfile(native_file), delete(native_file); end
    if isfolder(work_dir), rmdir(work_dir); end
end


function paired = fit_paired_sensitivities(DCMs, subjects, groups, FD, ...
        discovery_index, Xwithin, PEBwithin, nodeIDs)
% Fit window-FD bridges and paired cohorts with the original mapped priors.
% Output: lower/upper PEBs and their cohort, design and prior metadata; no files.
    subjects = string(subjects(:));
    groups = double(groups(:));
    discovery_index = discovery_index(:);
    n_reference = numel(discovery_index);
    group_center = mean(groups(discovery_index));
    n_edges = numel(nodeIDs)^2;
    months = [0,3,8];
    Xpair = [1,-0.5; 1,0.5];
    Dwithin = diag(size(Xwithin,1)./sum(Xwithin.^2,1));
    available = ~cellfun(@isempty, DCMs);

    paired = struct;
    paired.Model_names = {'Three_visit_windowFD', 'Paired_complete', 'Paired_all_eligible'};
    paired.End_months = months(2:3);
    paired.Reference_subjects = subjects(discovery_index);
    paired.Group_center = group_center;
    paired.Node_IDs = string(nodeIDs(:));
    paired.X_within = Xwithin;
    paired.X_pair = Xpair;
    paired.Within_prior_scaling = Dwithin;
    paired.Lower = cell(2,1);
    paired.Upper = cell(2,3);
    paired.Cohorts = cell(2,1);
    paired.Reference_pair_rows = cell(2,1);
    paired.X_upper = cell(2,3);
    paired.Input_M_lower = cell(2,1);
    paired.Input_M_upper = cell(2,3);
    paired.Temporal_mapping = cell(2,1);
    paired.Parameter_mapping = cell(2,1);
    paired.Pair_prior_scaling = cell(2,1);
    paired.Mapped_lower_prior_mean = cell(2,1);
    paired.Mapped_lower_prior_covariance = cell(2,1);
    paired.FD_center = cell(2,1);
    paired.FD_scale = cell(2,1);
    paired.Upper_prior_scaling = cell(2,1);

    for w = 1:2
        visit = w+1;
        visits = [1,visit];
        valid_motion = all(isfinite(FD(:,visits)) & FD(:,visits) >= 0, 2);
        rows = find(all(available(:,visits),2) & valid_motion);
        [~, reference_rows] = ismember(discovery_index, rows);
        paired.Reference_pair_rows{w} = reference_rows;

        % Keep the discovery motion scale when enlarging to all available pairs.
        motion_reference = FD(discovery_index,visits);
        fd_center = mean(motion_reference,1);
        fd_scale = std(motion_reference,1,1);
        fd_z = (FD(rows,visits)-fd_center)./fd_scale;
        Xreference = [ones(n_reference,1), groups(discovery_index)-group_center, fd_z(reference_rows,:)];
        Xall = [ones(numel(rows),1), groups(rows)-group_center, fd_z];
        bX = diag(n_reference./sum(Xreference.^2,1));
        paired.FD_center{w} = fd_center;
        paired.FD_scale{w} = fd_scale;
        paired.Upper_prior_scaling{w} = bX;
        paired.Cohorts{w} = table(rows, subjects(rows), groups(rows), ...
            ismember(rows,discovery_index), FD(rows,1), FD(rows,visit), ...
            fd_z(:,1), fd_z(:,2), 'VariableNames', ...
            {'Input_row','Participant_ID','Group_code','In_reference', ...
             'FD0','FDend','FD0_z','FDend_z'});

        % Match the pair mean/change prior to the original three-visit basis.
        T = [mean(Xwithin(visits,:),1); Xwithin(visit,:)-Xwithin(1,:)];
        K = kron(T, eye(n_edges));
        Dpair = T*Dwithin*T';
        paired.Temporal_mapping{w} = T;
        paired.Parameter_mapping{w} = K;
        paired.Pair_prior_scaling{w} = Dpair;
        paired.Mapped_lower_prior_mean{w} = K*spm_vec(PEBwithin{1}.M.pE);
        paired.Mapped_lower_prior_covariance{w} = K*PEBwithin{1}.M.pC*K';
        M = struct('X',Xpair, 'Xnames',{{'Pair_mean','Endpoint_minus_baseline'}}, ...
            'bX',Dpair, 'alpha',1, 'beta',16, 'Q','single', ...
            'hE',0, 'hC',1/16, 'maxit',256);
        paired.Input_M_lower{w} = M;
        paired.Lower{w} = cell(numel(rows),1);
        for i = 1:numel(rows)
            pair = DCMs(rows(i),visits)';
            paired.Lower{w}{i} = spm_dcm_peb(pair, M, {'A'});
        end

        % Bridge, complete-pair and enlarged-pair fits share reference scaling.
        Psets = {PEBwithin(:), paired.Lower{w}(reference_rows), paired.Lower{w}};
        Xsets = {Xreference, Xreference, Xall};
        for m = 1:3
            if m == 1
                first = PEBwithin{1};
            else
                first = paired.Lower{w}{reference_rows(1)};
            end
            pe = spm_vec(first.M.pE);
            pc = full(first.M.pC);
            M = struct('X',Xsets{m}, 'Xnames',{{'Mean','Group','FD0','FDend'}}, ...
                'bE',pe, 'bC',pc, 'pC',pc/16, 'bX',bX, ...
                'Q','single', 'hE',0, 'hC',1/16, 'maxit',256);
            paired.X_upper{w,m} = Xsets{m};
            paired.Input_M_upper{w,m} = M;
            paired.Upper{w,m} = spm_dcm_peb(Psets{m}, M, (1:numel(pe))');
        end
    end
end


function [coefficients,contrasts,edge_dictionary,scan_table] = make_peb_tables(PEB)
% Aim: prepare scientific tables from completed PEBs and original DCM summaries.
% Method: full-Cp contrasts and separate BMA estimates; output: four tables, no files.

ids = string(PEB.Node_IDs(:));
labels = string(PEB.Node_labels(:));
n_nodes = numel(ids); n_A = n_nodes^2;
n_within = size(PEB.Xwithin,2);
[target,source] = ind2sub([n_nodes,n_nodes],(1:n_A)');
extrinsic = find(target ~= source);
csv_columns = strings(n_A,1);
for p = 1:n_A
    csv_columns(p) = string(matlab.lang.makeValidName(sprintf('Edge_%s_to_%s_Hz', ...
        char(ids(source(p))),char(ids(target(p))))));
end

main_models = {PEB.NoFD,PEB.FullFD};
main_bma = {PEB.BMA_NoFD,PEB.BMA_FD};
main_names = ["No_FD","FD_adjusted"];
bma_mean = cell(1,2); bma_pp = cell(1,2);
selected = false(n_A,2);
group_terms = n_within + [2,3];
for m = 1:2
    bma_mean{m} = reshape(full(main_bma{m}.Ep),n_A,[]);
    bma_pp{m} = reshape(full(main_bma{m}.Pp),n_A,[]);
    selected(extrinsic,m) = any(bma_pp{m}(extrinsic,group_terms) > PEB.Threshold,2);
end
selected_union = find(any(selected,2));
discovery_groups = PEB.Groups(PEB.Discovery_index);
main_counts = [numel(discovery_groups),sum(discovery_groups == -1),sum(discovery_groups == 1)];
term_names = ["Group_x_Linear","Group_x_Quadratic"];
term_units = ["Hz/month","Hz/month^2"];

% Report both treatment-time coefficients for every extrinsic edge in both main models.
coefficient_rows = cell(numel(extrinsic)*2*2,24); r = 0;
for m = 1:2
    model = main_models{m};
    ep = model.Ep(:);
    for term = 1:2
        coefficient_column = group_terms(term);
        for p = extrinsic'
            index = (coefficient_column-1)*n_A+p;
            native_mean = full(ep(index));
            estimate = 2*native_mean;
            sd = 2*sqrt(full(model.Cp(index,index)));
            r = r+1;
            coefficient_rows(r,:) = {main_names(m),p,ids(source(p)),ids(target(p)), ...
                labels(source(p)),labels(target(p)),csv_columns(p),term_names(term), ...
                term_units(term),native_mean,2,estimate,estimate-1.959963984540054*sd, ...
                estimate+1.959963984540054*sd,2*bma_mean{m}(p,coefficient_column), ...
                bma_pp{m}(p,coefficient_column),selected(p,1),selected(p,2), ...
                all(selected(p,:)),main_counts(1),main_counts(2),main_counts(3),PEB.Threshold, ...
                "Full-PEB pointwise 95% CrI; not BMA or selection-adjusted"};
        end
    end
end
coefficients = cell2table(coefficient_rows,'VariableNames', ...
    {'Model','A_Linear_Index','Source_Node','Target_Node','Source_label','Target_label', ...
     'CSV_Column','Term','Unit','Native_fullPEB_half_difference','Group_scale_factor', ...
     'FullPEB_mean','FullPEB_CrI95_low','FullPEB_CrI95_high','BMA_MeditationMinusControl', ...
     'BMA_Pp','Selected_No_FD','Selected_FD_adjusted','Selected_both', ...
     'N','Control_N','Meditation_N','Pp_threshold','Interval'});

% Retain stable EC identifiers; corrected BMA column names explicitly denote full differences.
dictionary_rows = cell(numel(selected_union),19);
for e = 1:numel(selected_union)
    p = selected_union(e);
    dictionary_rows(e,:) = {p,ids(source(p)),ids(target(p)),csv_columns(p), ...
        labels(source(p)),labels(target(p)),selected(p,1),selected(p,2),all(selected(p,:)), ...
        2*bma_mean{1}(p,group_terms(1)),bma_pp{1}(p,group_terms(1)), ...
        2*bma_mean{1}(p,group_terms(2)),bma_pp{1}(p,group_terms(2)), ...
        2*bma_mean{2}(p,group_terms(1)),bma_pp{2}(p,group_terms(1)), ...
        2*bma_mean{2}(p,group_terms(2)),bma_pp{2}(p,group_terms(2)), ...
        PEB.Threshold,"BMA inclusion support within the respective reduction procedure"};
end
edge_dictionary = cell2table(dictionary_rows,'VariableNames', ...
    {'A_Linear_Index','Source_Node','Target_Node','CSV_Column','Source_label','Target_label', ...
     'Selected_No_FD','Selected_FD_adjusted','Selected_both', ...
     'NoFD_linear_BMA_MeditationMinusControl_Hz_per_month','NoFD_linear_Pp', ...
     'NoFD_quadratic_BMA_MeditationMinusControl_Hz_per_month2','NoFD_quadratic_Pp', ...
     'FD_linear_BMA_MeditationMinusControl_Hz_per_month','FD_linear_Pp', ...
     'FD_quadratic_BMA_MeditationMinusControl_Hz_per_month2','FD_quadratic_Pp', ...
     'Pp_threshold','Pp_interpretation'});

% Five specifications per window use their own full posterior covariance.
contrast_rows = cell(numel(extrinsic)*(2*5+3*2*2),23); r = 0;
for w = 1:2
    visit = w+1; month = PEB.Time(visit);
    cohort_groups = PEB.Paired.Cohorts{w}.Group_code;
    reference_groups = cohort_groups(PEB.Paired.Reference_pair_rows{w});
    paired_counts = [numel(reference_groups),sum(reference_groups == -1),sum(reference_groups == 1)];
    all_counts = [numel(cohort_groups),sum(cohort_groups == -1),sum(cohort_groups == 1)];
    window_motion = "FD0, FD"+string(month)+"; scaled on discovery cohort";
    for specification = 1:5
        if specification <= 2
            model = main_models{specification}; bma = main_bma{specification};
            name = main_names(specification); model_id = name;
            delta_basis = PEB.Xwithin(visit,:)-PEB.Xwithin(1,:);
            counts = main_counts;
            if specification == 1
                motion = "None";
            else
                motion = "FD0, FD3, FD8; scaled on discovery cohort";
            end
        else
            pair_model = specification-2;
            model = PEB.Paired.Upper{w,pair_model}; bma = [];
            name = string(PEB.Paired.Model_names{pair_model});
            model_id = name+"_0_"+string(month);
            motion = window_motion;
            if pair_model == 1
                delta_basis = PEB.Xwithin(visit,:)-PEB.Xwithin(1,:);
                counts = main_counts;
            else
                delta_basis = PEB.Paired.X_pair(2,:)-PEB.Paired.X_pair(1,:);
                if pair_model == 2, counts = paired_counts; else, counts = all_counts; end
            end
        end
        n_upper = numel(model.Ep)/(n_A*numel(delta_basis));
        upper_row = zeros(1,n_upper); upper_row(2) = 2;
        for p = extrinsic'
            weight = peb_table_weight(p,delta_basis,upper_row,n_A);
            r = r+1;
            contrast_rows(r,:) = peb_table_contrast(model,bma,weight,name,model_id, ...
                "Group_difference_in_change","meditation-minus-control",month,p, ...
                source,target,ids,labels,selected,counts,motion);
        end
    end
end

% Fitted main-model trajectories combine average, linear and quadratic components.
% Motion-adjusted means set all standardized FD covariates to zero for both arms.
group_names = ["control","meditation"]; group_values = [-1,1];
for m = 1:2
    model = main_models{m}; bma = main_bma{m};
    n_upper = numel(model.Ep)/(n_A*n_within);
    if m == 1
        motion = "None";
    else
        motion = "FD0=FD3=FD8=0 on the discovery-standardized scale";
    end
    for visit = 1:numel(PEB.Time)
        for g = 1:2
            upper_row = zeros(1,n_upper);
            upper_row(1:2) = [1,group_values(g)-PEB.Group_center];
            for p = extrinsic'
                weight = peb_table_weight(p,PEB.Xwithin(visit,:),upper_row,n_A);
                r = r+1;
                contrast_rows(r,:) = peb_table_contrast(model,bma,weight,main_names(m),main_names(m), ...
                    "Fitted_group_mean",group_names(g),PEB.Time(visit),p, ...
                    source,target,ids,labels,selected,main_counts,motion);
            end
        end
    end
end
contrasts = cell2table(contrast_rows,'VariableNames', ...
    {'Model','Model_ID','Contrast','Group','End_month','A_Linear_Index', ...
     'Source_Node','Target_Node','Source_label','Target_label','FullPEB_mean_Hz', ...
     'FullPEB_CrI95_low_Hz','FullPEB_CrI95_high_Hz','BMA_mean_Hz', ...
     'Selected_No_FD','Selected_FD_adjusted','Selected_both','N','Control_N', ...
     'Meditation_N','Motion_covariates','Unit','Interval'});

% Export original scan means without clinical or trajectory-label exclusions.
subjects = string(PEB.Subjects(:)); groups = PEB.Groups(:);
in_discovery = false(numel(subjects),1); in_discovery(PEB.Discovery_index) = true;
scan_rows = cell(0,5+numel(selected_union));
for i = 1:numel(subjects)
    group = group_names(1+(groups(i) == 1));
    for visit = 1:numel(PEB.Time)
        dcm = PEB.DCM_summaries{i,visit};
        if isempty(dcm), continue; end
        original_A = full(dcm.Ep.A);
        original_values = reshape(original_A(selected_union),1,[]);
        scan_rows(end+1,:) = [{subjects(i),PEB.Time(visit),group,in_discovery(i),PEB.FD(i,visit)}, ...
            num2cell(original_values)];
    end
end
scan_table = cell2table(scan_rows,'VariableNames', ...
    [{'Participant_ID','Session_Month','Group','In_PEB_Discovery','mean_fd'}, ...
     cellstr(csv_columns(selected_union))']);
end

function weight = peb_table_weight(edge,within_row,upper_row,n_A)
% SPM vector order is edge within time component within upper covariate.
values = kron(upper_row,within_row);
weight = zeros(n_A*numel(values),1);
weight(edge:n_A:end) = values(:);
end

function row = peb_table_contrast(model,bma,weight,name,model_id,contrast,group,month, ...
        edge,source,target,ids,labels,selected,counts,motion)
% Full-PEB means and uncertainty share the same model; BMA has no exported interval.
estimate = full(weight'*model.Ep(:));
sd = sqrt(full(weight'*model.Cp*weight));
bma_estimate = NaN;
if ~isempty(bma), bma_estimate = full(weight'*bma.Ep(:)); end
row = {name,model_id,contrast,group,month,edge,ids(source(edge)),ids(target(edge)), ...
    labels(source(edge)),labels(target(edge)),estimate, ...
    estimate-1.959963984540054*sd,estimate+1.959963984540054*sd,bma_estimate, ...
    selected(edge,1),selected(edge,2),all(selected(edge,:)),counts(1),counts(2),counts(3), ...
    motion,"Hz","Full-PEB pointwise 95% CrI; not BMA or selection-adjusted"};
end
