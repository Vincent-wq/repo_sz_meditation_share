function peb_section4_final
% Estimate clinical-neural correspondence from the original fitted DCMs.
% Fit the eight Section 3 edges, with two primary edges identified for reporting.
% Use the formal Section 3 edge dictionary and the Python clinical sample.
% Save clinical contrasts, sensitivity summaries and temporal-block evidence as CSVs.

%% 1. Files and model settings
project_dir = 'C:\PRIVATE_INPUTS\repo_sz_meditation';
dcm_dir = ['C:\PRIVATE_INPUTS\sz_meditation\derivatives\' ...
    'results_20260701\spDCM_Estimated_TR1p5_TE30'];
output_dir = fullfile(project_dir, 'data', 'neural_clinical');
clinical_file = fullfile(output_dir, 'Neural_clinical_sample.csv');
edge_file = fullfile(project_dir, 'data', 'dcm', 'PEB_section3', ...
    'Selected_edge_dictionary.csv');
run_omissions = true; % Include all participant omissions for the influence table.

spm('defaults', 'FMRI');
spm_get_defaults('dcm.verbose', true);
if ~isfolder(output_dir), mkdir(output_dir); end
clinical = readtable(clinical_file, 'TextType', 'string');
edges = sortrows(readtable(edge_file, 'TextType', 'string'), 'A_Linear_Index');
covariates = {'PANSS_SUM_base','RBANS_sum_base','CPZ_dose','FD0','FD3','FD8'};
outcomes = {'PANSS_dec_rate_3','RBANS_inc_rate_3', ...
    'PANSS_dec_rate_8','RBANS_inc_rate_8'};
required = [{'participant_id','treatment','PANSS_SUM_3','RBANS_sum_3'}, ...
    covariates,outcomes];
if ~all(ismember(required, clinical.Properties.VariableNames))
    error('PEB:Input', 'Run the Python input block first; required clinical columns are missing.');
end
ids = strtrim(string(clinical.participant_id));
n = height(clinical);
values = clinical{:, [covariates,outcomes]};
raw_scores = clinical{:, {'PANSS_SUM_3','RBANS_sum_3'}};
if numel(unique(ids)) ~= n || ~all(ismember(clinical.treatment,[0,1]))
    error('PEB:Sample', 'The fresh sample must contain one row per participant and treatment coded 0/1.');
end
if any(~isfinite(values(:))) || any(~isfinite(raw_scores(:))) || ...
        any(std(values,1,1) == 0) || any(std(raw_scores,1,1) == 0)
    error('PEB:Input', 'Clinical and FD fields must be finite and nonconstant; no imputation is applied.');
end

% Match the Python population-SD convention on this same complete-case cohort.
z = (values-mean(values,1))./std(values,1,1);
raw_z = (raw_scores-mean(raw_scores,1))./std(raw_scores,1,1);
raw_z(:,1) = -raw_z(:,1); % Higher values indicate lower PANSS3 at fixed baseline.
group = 2*clinical.treatment-1;
core = [ones(n,1),group-mean(group),z(:,1:6)];
core_names = [{'Mean','Group'},covariates];

%% 2. The eight Section 3 edges and three-visit neural representation
% Select the declared eight before fitting; no clinical-outcome selection is used.
% The complete original DCM network remains the source of first-stage uncertainty.
original_eight = ismember(lower(strtrim(string(edges.Selected_No_FD))), ...
    ["1","1.0","true"]);
edges = edges(original_eight,:);
node_ids = {'IC00_SMN';'IC02_SalVentA';'IC04_Vis_V1';'IC05_aDMN'; ...
    'IC09_pDMN';'IC15_DAN';'IC16_FPCN';'IC18_SalVentB';'IC19_ECN'};
a = edges.A_Linear_Index(:);
if height(edges) ~= 8 || numel(unique(a)) ~= 8 || any(a < 1 | a > 81)
    error('PEB:Edges', 'The correspondence specification requires the eight Section 3 edges marked Selected_No_FD.');
end
[target,source] = ind2sub([9,9],a);
if any(target == source) || ...
        ~isequal(string(edges.Source_Node),string(node_ids(source))) || ...
        ~isequal(string(edges.Target_Node),string(node_ids(target)))
    error('PEB:Edges', 'The edge dictionary does not match the established A(target,source) node order.');
end
edges.Original_eight = true(height(edges),1);
edges.Primary_edge = ismember(a,[10,17]); % IC02 -> IC00 and IC02 -> IC18.
q = [a;a+81;a+162];
n_parameters = numel(q); % Eight edges x three within-person time coefficients.
months = [0;3;8];
linear_time = months-mean(months);
quadratic_time = linear_time.^2-linear_time* ...
    ((linear_time'*linear_time)\(linear_time'*(linear_time.^2)));
quadratic_time = quadratic_time-mean(quadratic_time);
Xwithin = [ones(3,1),linear_time,quadratic_time];

files = cell(n,3);
for i = 1:n
    for j = 1:3
        files{i,j} = fullfile(dcm_dir,sprintf('DCM_%s_ses-%d_timecourses.mat', ...
            char(ids(i)),months(j)));
        if ~isfile(files{i,j})
            error('PEB:DCM', 'Required estimated DCM is missing: %s',files{i,j});
        end
    end
end

% Native SPM may write tmp.mat internally; keep it outside the result directory.
original_dir = pwd;
work_dir = tempname;
mkdir(work_dir);
restore_dir = onCleanup(@() finish_work(original_dir,work_dir));
cd(work_dir);
lower_peb = cell(n,1);
for i = 1:n
    lower_peb{i} = fit_peb(files(i,:)',Xwithin,{'A'});
    if ~isequal(size(lower_peb{i}.Ep),[81,3]) || ...
            ~isequal(lower_peb{i}.Pind(:),(1:81)')
        error('PEB:Parameters', 'The fitted DCM A-parameter order differs for %s.',char(ids(i)));
    end
    fprintf('Within-person PEB %s: %d/%d\n',char(ids(i)),i,n);
end
fprintf('Matched N=%d: %d control/%d meditation; original DCM estimates retained.\n', ...
    n,sum(clinical.treatment == 0),sum(clinical.treatment == 1));

%% 3. Joint, conditioning and clinical-definition models
% Each clinical coefficient predicts neural parameters;
designs = {[core,z(:,7:8)], [core,z(:,7:10)], ...
    [core,z(:,7)], [core,z(:,8)], [core,raw_z]};
names = {[core_names,outcomes(1:2)], [core_names,outcomes], ...
    [core_names,outcomes(1)], [core_names,outcomes(2)], ...
    [core_names,{'PANSS3_lower_z','RBANS3_higher_z'}]};
model_names = {'Early_correspondence','Late_conditional_correspondence', ...
    'PDR3_only','RIR3_only','Raw_score_early'};
models = cell(5,1);
for m = 1:5
    M = upper_design(designs{m},names{m});
    models{m} = fit_peb(lower_peb,M,q);
    if ~isequal(models{m}.Pind(:),q) || ...
            ~isequal(size(models{m}.Ep),[n_parameters,size(designs{m},2)])
        error('PEB:Parameters', 'Unexpected selected-parameter order in %s.',model_names{m});
    end
    fprintf('%s: completed.\n',model_names{m});
end

%% 4. Full-model clinical contrasts and the two sensitivity definitions
% Transform the whole posterior covariance to 0-3 or 0-8 EC-change contrasts.
early = contrast_table(models{1},names{1},outcomes(1:2),Xwithin(2,:)-Xwithin(1,:), ...
    model_names{1},3,edges,n);
late = contrast_table(models{2},names{2},outcomes(3:4),Xwithin(3,:)-Xwithin(1,:), ...
    model_names{2},8,edges,n);
correspondence = [early;late];
correspondence = arm_counts(correspondence,clinical.treatment);
writetable(correspondence,fullfile(output_dir,'PEB_clinical_correspondence.csv'));

symptom_only = contrast_table(models{3},names{3},outcomes(1),Xwithin(2,:)-Xwithin(1,:), ...
    model_names{3},3,edges,n);
cognition_only = contrast_table(models{4},names{4},outcomes(2),Xwithin(2,:)-Xwithin(1,:), ...
    model_names{4},3,edges,n);
conditioning = [early;symptom_only;cognition_only];
conditioning.Clinical_conditioning = repmat("Mutually adjusted PDR3 and RIR3",height(conditioning),1);
conditioning.Clinical_conditioning(conditioning.Model ~= "Early_correspondence") = "One clinical domain";
conditioning = arm_counts(conditioning,clinical.treatment);
writetable(conditioning,fullfile(output_dir,'PEB_clinical_conditioning.csv'));

raw = contrast_table(models{5},names{5},names{5}(9:10),Xwithin(2,:)-Xwithin(1,:), ...
    model_names{5},3,edges,n);
definition = [early;raw];
definition.Clinical_definition = [repmat("Ratio",height(early),1);repmat("Raw_score",height(raw),1)];
definition = arm_counts(definition,clinical.treatment);
writetable(definition,fullfile(output_dir,'PEB_clinical_raw_score.csv'));

%% 5. Evidence for the two complete early clinical-temporal blocks
% Use actual full-model priors and Gaussian posterior densities at zero.
% Visit-average clinical effects and all nuisance effects remain free.
temporal = temporal_table(models([1,5]),names([1,5]),n,edges);
temporal = arm_counts(temporal,clinical.treatment);
writetable(temporal,fullfile(output_dir,'PEB_clinical_temporal_blocks.csv'));

%% 6. Participant influence on the early joint correspondence
% Keep the full-sample scales/design entries; omit one participant without reselection.
% The inherited native prior rule remains fixed, not its numerical N-dependent values.
if run_omissions
    values_loo = nan(n,height(early),3);
    omission_status = repmat("Failed",n,1);
    omission_error = repmat("",n,1);
    for i = 1:n
        keep = (1:n)' ~= i;
        try
            M = upper_design(designs{1}(keep,:),names{1});
            fit = fit_peb(lower_peb(keep),M,q);
            row = contrast_table(fit,names{1},outcomes(1:2),Xwithin(2,:)-Xwithin(1,:), ...
                model_names{1},3,edges,n-1);
            values_loo(i,:,:) = reshape(row{:, ...
                {'FullPEB_Hz_per_clinical_SD','CrI95_low','CrI95_high'}},1,height(early),3);
            omission_status(i) = "Usable"; % Native fit and contrast calculation completed.
        catch problem
            omission_error(i) = string(problem.message);
        end
        fprintf('Omit %s: %s (%d/%d)\n',char(ids(i)),char(omission_status(i)),i,n);
    end
    influence = influence_table(early,values_loo,ids,omission_status,omission_error);
else
    influence = influence_table(early,nan(n,height(early),3),ids,repmat("Not_run",n,1), ...
        repmat("Omission module disabled",n,1));
    fprintf('Participant omissions were not run; the influence table marks them Not_run.\n');
end
influence = arm_counts(influence,clinical.treatment);
for i = 1:n
    row = influence.Omitted_ID == ids(i);
    influence.Control_N(row) = influence.Control_N(row)-(clinical.treatment(i) == 0);
    influence.Meditation_N(row) = influence.Meditation_N(row)-(clinical.treatment(i) == 1);
end
writetable(influence,fullfile(output_dir,'PEB_clinical_LOO.csv'));
clear restore_dir
disp(correspondence);
disp(temporal);
fprintf('Results CSVs saved in: %s\n',output_dir);
fprintf('Eight ECs were fitted jointly; the two primary edges are reporting labels within this model.\n');
fprintf('All three neural visits and FD0/3/8 inform the retrospective correspondence models.\n');
fprintf('Clinical slopes are Hz per clinical SD; do not apply the group-contrast factor of two.\n');
fprintf('Intervals are full-model pointwise CrIs; omission ranges are influence summaries, not CIs.\n');
end

%% functions
function M = upper_design(X,names)
% Keep the original inherited PEB prior and precision settings.
    if rank(X) < size(X,2)
        error('PEB:Design','The clinical design is rank deficient.');
    end
    M = struct('X',X,'Xnames',{names},'maxit',256, ...
        'alpha',1,'beta',16,'Q','single','hE',0,'hC',1/16);
end


function fit = fit_peb(P,M,parameters)
% Fit with native SPM settings and return finite posterior estimates.
    fit = spm_dcm_peb(P,M,parameters);
    if any(~isfinite(spm_vec(fit.Ep))) || any(~isfinite(fit.Cp(:))) || ~isfinite(fit.F)
        error('PEB:Fit','SPM returned a nonfinite posterior or free energy.');
    end
end


function table_out = contrast_table(fit,names,domains,w,label,month,edges,n)
% One row per selected edge and clinical regressor, retaining all pointwise intervals.
    ne = height(edges);
    np = 3*ne;
    rows = cell(ne*numel(domains),16);
    r = 0;
    for c = 1:numel(domains)
        column = find(strcmp(names,domains{c}));
        for e = 1:ne
            ix = (column-1)*np+[e,e+ne,e+2*ne];
            estimate = full(w*reshape(fit.Ep(ix),3,1));
            variance = full(w*fit.Cp(ix,ix)*w');
            if ~isfinite(variance) || variance <= 0
                error('PEB:Contrast','A requested contrast has nonpositive or invalid variance.');
            end
            halfwidth = 1.959963984540054*sqrt(variance);
            r = r+1;
            rows(r,:) = {string(label),string(domains{c}),"0-"+string(month)+" months",month, ...
                edges.A_Linear_Index(e),string(edges.Source_Node(e)),string(edges.Target_Node(e)), ...
                string(edges.Source_label(e)),string(edges.Target_label(e)),edges.Original_eight(e), ...
                estimate,estimate-halfwidth,estimate+halfwidth,n, ...
                "Hz EC change per one full-sample clinical SD", ...
                "Full-PEB pointwise 95% CrI; not selection-adjusted"};
        end
    end
    table_out = cell2table(rows,'VariableNames',{'Model','Clinical_regressor', ...
        'EC_change_window','End_month','A_Linear_Index','Source_Node','Target_Node', ...
        'Source_label','Target_label','Original_eight','FullPEB_Hz_per_clinical_SD', ...
        'CrI95_low','CrI95_high','N','Unit','Interval'});
    table_out.Primary_edge = repmat(edges.Primary_edge,numel(domains),1);
    table_out.EC_count = repmat(ne,height(table_out),1);
    table_out.Edge_set = repmat("Section3_original_eight",height(table_out),1);
end


function table_out = temporal_table(models,design_names,n,edges)
% Exact Gaussian zero-constraint evidence within each four-model clinical family.
    patterns = logical([0,0;1,0;0,1;1,1]);
    labels = ["Neither_temporal","Symptoms_temporal","Cognition_temporal","Both_temporal"];
    definitions = ["Ratio","Raw_score"];
    ne = height(edges);
    np = 3*ne;
    clinical_names = {{'PANSS_dec_rate_3','RBANS_inc_rate_3'}, ...
        {'PANSS3_lower_z','RBANS3_higher_z'}};
    rows = cell(8,12);
    r = 0;
    for d = 1:2
        fit = models{d};
        qe = spm_vec(fit.Ep); qc = full(fit.Cp);
        pe = spm_vec(fit.M.pE); pc = full(fit.M.pC);
        full_dimension = numel(qe);
        if numel(pe) ~= full_dimension || ~isequal(size(pc),[full_dimension,full_dimension])
            error('PEB:Prior','Temporal evidence requires the actual complete prior in the fitted PEB.');
        end
        block = cell(1,2);
        for c = 1:2
            column = find(strcmp(design_names{d},clinical_names{d}{c}));
            block{c} = (column-1)*np+(ne+1:np);
        end
        delta = zeros(1,4);
        for k = 1:4
            dropped = [block{~patterns(k,:)}];
            if isempty(dropped), continue; end
            retained = setdiff(1:full_dimension,dropped);
            cross = pc(dropped,retained);
            scale = sqrt(norm(pc(dropped,dropped),'fro')*norm(pc(retained,retained),'fro'));
            if norm(cross,'fro') > 1e-10*max(eps,scale) || any(abs(pe(dropped)) > 1e-12)
                error('PEB:Prior','The intended zero blocks are coupled to retained parameters in the full prior.');
            end
            Rq = chol((qc(dropped,dropped)+qc(dropped,dropped)')/2);
            Rp = chol((pc(dropped,dropped)+pc(dropped,dropped)')/2);
            zq = Rq'\qe(dropped);
            zp = Rp'\pe(dropped);
            delta(k) = sum(log(diag(Rp)))-sum(log(diag(Rq)))+0.5*(zp'*zp-zq'*zq);
        end
        if any(~isfinite(delta)), error('PEB:Evidence','Nonfinite temporal-block log evidence.'); end
        weight = exp(delta-max(delta));
        probability = weight/sum(weight);
        inclusion = zeros(1,2); inclusion_logBF = zeros(1,2);
        for c = 1:2
            inclusion(c) = sum(probability(patterns(:,c)));
            inclusion_logBF(c) = logsumexp(delta(patterns(:,c)))-logsumexp(delta(~patterns(:,c)));
        end
        for k = 1:4
            r = r+1;
            rows(r,:) = {definitions(d),labels(k),n,patterns(k,1),patterns(k,2), ...
                0.25,delta(k),probability(k),inclusion(1),inclusion(2), ...
                inclusion_logBF(1),inclusion_logBF(2)};
        end
    end
    table_out = cell2table(rows,'VariableNames',{'Clinical_definition','Model','N', ...
        'Symptom_temporal_block','Cognition_temporal_block','Prior_model_probability', ...
        'LogBF_vs_both_temporal','Posterior_model_probability', ...
        'Symptom_inclusion_probability','Cognition_inclusion_probability', ...
        'Symptom_inclusion_logBF','Cognition_inclusion_logBF'});
    table_out.Method = repmat("Gaussian posterior/prior density ratio at zero; within-definition comparison", ...
        height(table_out),1);
    table_out.EC_count = repmat(ne,height(table_out),1);
    table_out.Temporal_coefficients_per_clinical_block = repmat(2*ne,height(table_out),1);
    table_out.Edge_set = repmat("Section3_original_eight",height(table_out),1);
end


function value = logsumexp(x)
    peak = max(x);
    value = peak+log(sum(exp(x-peak)));
end


function table_out = influence_table(reference,values,ids,status,messages)
% Reference summaries and every omission stay together; none is silently discarded.
    n = numel(ids); nr = height(reference);
    table_out = repmat(reference,n+1,1);
    count = height(table_out);
    table_out.Omitted_ID = repmat("",count,1);
    table_out.Fit_status = repmat("Reference",count,1);
    table_out.Error = repmat("",count,1);
    table_out.Reference_estimate = repmat(reference.FullPEB_Hz_per_clinical_SD,n+1,1);
    table_out.LOO_usable_min = nan(count,1);
    table_out.LOO_usable_max = nan(count,1);
    table_out.LOO_usable_opposite_sign_count = nan(count,1);
    usable = status == "Usable";
    table_out.Usable_omissions = repmat(sum(usable),count,1);
    table_out.Failed_omissions = repmat(sum(status == "Failed"),count,1);
    table_out.Not_run_omissions = repmat(sum(status == "Not_run"),count,1);
    table_out.Total_planned_omissions = repmat(n,count,1);
    table_out.All_omissions_usable = repmat(all(usable),count,1);
    for r = 1:nr
        v = reshape(values(usable,r,1),[],1);
        if ~isempty(v)
            table_out.LOO_usable_min(r) = min(v);
            table_out.LOO_usable_max(r) = max(v);
            table_out.LOO_usable_opposite_sign_count(r) = sum(v*reference.FullPEB_Hz_per_clinical_SD(r) < 0);
        end
    end
    for i = 1:n
        rows = i*nr+(1:nr);
        table_out{rows,{'FullPEB_Hz_per_clinical_SD','CrI95_low','CrI95_high'}} = reshape(values(i,:,:),nr,3);
        table_out.N(rows) = n-1;
        table_out.Omitted_ID(rows) = ids(i);
        table_out.Fit_status(rows) = status(i);
        table_out.Error(rows) = messages(i);
    end
    table_out.LOO_interpretation = repmat("Influence ranges among usable omissions, not confidence intervals",count,1);
end


function table_out = arm_counts(table_out,treatment)
% Keep the actual fresh cohort counts in each scientific result table.
    table_out.Control_N = repmat(sum(treatment == 0),height(table_out),1);
    table_out.Meditation_N = repmat(sum(treatment == 1),height(table_out),1);
end


function finish_work(original_dir,work_dir)
% Restore the directory; remove only this run's native temporary file.
    cd(original_dir);
    filename = fullfile(work_dir,'tmp.mat');
    if isfile(filename), delete(filename); end
    if isfolder(work_dir), rmdir(work_dir); end
end
