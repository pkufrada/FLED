% EFIRM_L858R_NSCLC_all.m
% Function: Perform fragment decomposition (Short/Medium/Long) for NSCLC EGFR L858R samples 
% and compare with tissue ddPCR results.
% Description: Input CSV contains only "mean" values per probe without error bars. 
% This script uses equal weighting and half-LOD processing.
% Output: A complete set of charts/CSVs consistent with the KRAS version + iterative R^2 dynamic GIF.

%% I/O
in_csv  = '20250917 All NSCLC Clinical .csv';   % Use current directory or provide an absolute path
timestamp = datestr(now,'yyyymmdd_HHMMSS');
outdir = fullfile(pwd, ['results_L858R_' timestamp]);
if ~exist(outdir,'dir'), mkdir(outdir); end

% === Automatic script backup to txt (Consistent with KRAS version) ===
currentScriptName = mfilename('fullpath');
[~, scriptBaseName] = fileparts(currentScriptName);
backupName = fullfile(outdir, ['log_' scriptBaseName, '.txt']);
try
    copyfile([currentScriptName, '.m'], backupName);
    fprintf('✔ Current script backed up to: %s\n', backupName);
catch
    warning('⚠ Unable to backup current script; mfilename may not recognize the script.');
end

%% Read Table & Automatic Column Recognition
T = readtable(in_csv);

% ====== PATCH: Assemble Unique Sample Names <SampleID>_<Type>_<GIT> ======
% Automatic positioning of required columns
idCandidates   = {'sampleid','sample_id','id','case','patient','accession','sample'};
typeCandidates = {'type','sampletype','specimen','specimentype','matrix','source'};
gitCandidates  = {'git','git_status','gitstatus','gitresult','git_res','gitposneg','git_pn'};
vnamesLower = lower(T.Properties.VariableNames);

cID   = find(ismember(vnamesLower, idCandidates), 1);
cType = find(ismember(vnamesLower, typeCandidates), 1);
cGIT  = find(ismember(vnamesLower, gitCandidates), 1);

if isempty(cID) || isempty(cType) || isempty(cGIT)
    error(['Required columns for sample naming not found.\n' ...
        'Expected at least one of the following:\n' ...
        '  SampleID: %s\n  Type: %s\n  GIT: %s\n(Case-insensitive, synonyms supported)'], ...
        strjoin(idCandidates, ', '), strjoin(typeCandidates, ', '), strjoin(gitCandidates, ', '));
end

rawID   = string(T{:,cID});
rawType = string(T{:,cType});
rawGIT  = string(T{:,cGIT});

% Normalize Type → P/M/S
normType = strings(size(rawType));
for i=1:numel(rawType)
    t = lower(strtrim(rawType(i)));
    if contains(t,'plasma') || startsWith(t,'p')
        normType(i) = "P";
    elseif contains(t,'mpe') || contains(t,'pleural') || startsWith(t,'m')
        normType(i) = "M";
    elseif contains(t,'saliva') || startsWith(t,'s')
        normType(i) = "S";
    else
        % Default to uppercase first letter if unknown
        if strlength(t) >= 1
            normType(i) = upper(extractBetween(t,1,1));
        else
            normType(i) = "U";
        end
    end
end

% Normalize GIT → G+/G-
normGIT = strings(size(rawGIT));
for i=1:numel(rawGIT)
    g = lower(strtrim(rawGIT(i)));
    if contains(g,'+') || contains(g,'pos')
        normGIT(i) = "G+";
    elseif contains(g,'-') || contains(g,'neg')
        normGIT(i) = "G-";
    elseif g == "g+" || g == "g-"
        normGIT(i) = upper(g);
    else
        % Default to "G?" if unknown
        normGIT(i) = "G?";
    end
end

sampleNames = strtrim(rawID) + "_" + normType + "_" + normGIT;

% De-duplication (Add _2/_3 suffix for duplicates; first occurrence remains unchanged)
sampleNames = make_unique_names(sampleNames);
% ====== END PATCH ======

% Probe Columns (Short/Medium/Long)
cShort  = find(ismember(vnamesLower, lower({'Short','Close','EGFR_Short','KRAS_Close_Avg'})), 1);
cMed    = find(ismember(vnamesLower, lower({'Medium','Moderate','EGFR_Medium','KRAS_Moderate_Avg'})), 1);
cLong   = find(ismember(vnamesLower, lower({'Long','Far','EGFR_Long','KRAS_Far_Avg'})), 1);
if any(cellfun(@isempty, {cShort,cMed,cLong}))
    error(['Failed to identify Short/Medium/Long columns. Ensure CSV contains naming like:\n' ...
        '  Short/Medium/Long OR Close/Moderate/Far OR EGFR_Short/EGFR_Medium/EGFR_Long']);
end

% ddPCR Columns (Optional)
ddCols = {'ddPCR','ddpcr','tissue_ddpcr','Tissue_ddPCR','KRAS_total_ddPCR','ddPCR_cpsuL','ddPCR_cps_uL'};
cDD = find(ismember(vnamesLower, lower(ddCols)), 1);

% Extract Values (cps/µL)
Short_avg  = T{:,cShort};
Med_avg    = T{:,cMed};
Long_avg   = T{:,cLong};
if ~isempty(cDD)
    ddPCR_vec = T{:,cDD};
else
    ddPCR_vec = nan(height(T),1);
end

%% —— LOD + Equal Weights + log1p Fitting Preparation ——
X_obs = [Short_avg(:), Med_avg(:), Long_avg(:)];   % N×3 matrix
nz = X_obs(X_obs>0);
if isempty(nz), error('All probe values are 0; cannot perform fitting.'); end
lod = max(min(nz,[],'omitnan'), 1e-4);
X_obs(X_obs==0) = 0.5*lod;

W_mat = ones(size(X_obs));             % Equal weighting
X_log = log1p(X_obs);
EFIRM_total_sum3 = sum(X_obs,2);

%% —— Iterative Fitting: X_log ≈ b + Y_abs * A^T —— 
rng(0);
max_iter = 50; tol = 1e-8; prev_loss = Inf;
A_iter = eye(3);
b_iter = zeros(3,1);
Y_prop_init = repmat([0.6 0.05 0.35], numel(EFIRM_total_sum3), 1);
Y_abs = max(Y_prop_init .* EFIRM_total_sum3, 1e-3);
lambda_ridge = 1e-2; delta = 0.5;
W1 = diag(W_mat(:,1)); W2 = diag(W_mat(:,2)); W3 = diag(W_mat(:,3));

% Iterative R² dynamic plot + GIF setup
fig_iter = figure('Units','pixels','Position',[100, 100, 1000, 500]); hold on;
title('R^2 and Proportion Change per Iteration (log1p, WLS)');
xlabel('Iteration'); yyaxis left; ylabel('R^2'); ylim([0 1]);
r2_plot = plot(nan, nan, 'b-', 'DisplayName', 'R^2');
yyaxis right; ylabel('Mean Proportion'); ylim([0 1]);
p1 = plot(nan, nan, 'r-', 'DisplayName', 'Short');
p2 = plot(nan, nan, 'g-', 'DisplayName', 'Medium');
p3 = plot(nan, nan, 'm-', 'DisplayName', 'Long');
legend show;
gifPath = fullfile(outdir,'iterative_R2_convergence.gif'); makeGif = true;

% Robust initialization for b and A
for j = 1:3
    yj = X_log(:,j); Wj = W_mat(:,j);
    Z  = [ones(size(Y_abs,1),1), Y_abs];  t = yj;  w = Wj;
    beta = [0; A_iter(j,:)'];
    for ir=1:10
        r = t - Z*beta; s = median(abs(r))+1e-8; z = r/s;
        hub = ones(size(z)); idx = abs(z)>delta; hub(idx) = delta./abs(z(idx));
        W_eff = diag(w .* hub);
        beta = (Z'*W_eff*Z + lambda_ridge*eye(4)) \ (Z'*W_eff*t);
    end
    b_iter(j)   = beta(1); A_iter(j,:) = beta(2:4)';
end

loss_history = zeros(max_iter,1);
R2_history   = zeros(max_iter,1);
proportion_history = zeros(max_iter,3);
opts = optimoptions('fmincon','Display','off','Algorithm','sqp','MaxIterations',200);
epsilon = 1e-4;  Aeq = [1,1,1];  beq = 1;  lb = [epsilon,epsilon,epsilon];  ub = [1,1,1];

for iter=1:max_iter
    % === REPLACE: Update matrix A ===
        Y = Y_abs;
        Z = [ones(size(Y,1),1), Y];                 % Add intercept
        I4 = eye(4); I4(1,1) = 0;                   % Do not regularize intercept
        for j = 1:3
            wj = W_mat(:,j);
            Wj = diag(wj);
            beta = (Z'*Wj*Z + lambda_ridge*I4) \ (Z'*Wj*X_log(:,j));
            b_iter(j)   = beta(1);
            A_iter(j,:) = beta(2:4)';
        end

    % Fixed A and b; solve for proportions p per sample (simplex with damping + entropy regularization)
    Y_abs_new = zeros(size(Y_abs));
    gamma = 1e-3;      % Entropy regularization weight to avoid edge solutions
    alpha = 0.7;       % Damping step size

    for i=1:size(X_log,1)
        xi_log = X_log(i,:)';  
        wi     = W_mat(i,:)';  
        ktot   = EFIRM_total_sum3(i);
        
        % Objective function with entropy regularization
        fun = @(p) obj_huber_p(p, A_iter, b_iter, xi_log, wi, ktot, delta, gamma);
        
        % Initial normalized solution
        p0 = Y_abs(i,:); 
        p0 = p0./max(sum(p0),eps);
        
        % Optimization
        p_opt = fmincon(fun, p0, [],[], Aeq, beq, lb, ub, [], opts);
        
        % Damping update: prevent large jumps to extreme values
        p_old = Y_abs(i,:)./max(sum(Y_abs(i,:)),eps);
        p_new = (1-alpha)*p_old + alpha*p_opt(:)';
        Y_abs_new(i,:) = ktot * p_new;
    end

    % Evaluation and Dynamic Plots
    X_fit_log_it = b_iter' + Y_abs_new*A_iter';
    SSres = sum(((X_log - X_fit_log_it).^2).*W_mat, 'all');
    SStot = sum(((X_log - mean(X_log,'all')).^2).*W_mat, 'all');
    R2 = 1 - SSres/max(SStot,1e-8);
    loss_history(iter) = sqrt(SSres/numel(X_log));
    R2_history(iter)   = R2;
    Y_prop_it = Y_abs_new ./ max(sum(Y_abs_new,2), eps);
    proportion_history(iter,:) = mean(Y_prop_it,1,'omitnan');

    set(r2_plot, 'XData', 1:iter, 'YData', R2_history(1:iter));
    set(p1, 'XData', 1:iter, 'YData', proportion_history(1:iter,1));
    set(p2, 'XData', 1:iter, 'YData', proportion_history(1:iter,2));
    set(p3, 'XData', 1:iter, 'YData', proportion_history(1:iter,3));
    drawnow;
    
    if makeGif
        frame = getframe(fig_iter);
        [A_gif, map_gif] = rgb2ind(frame2im(frame), 256);
        if iter==1
            imwrite(A_gif, map_gif, gifPath, 'gif', 'LoopCount', Inf, 'DelayTime', 0.15);
        else
            imwrite(A_gif, map_gif, gifPath, 'gif', 'WriteMode', 'append', 'DelayTime', 0.15);
        end
    end

    if abs(prev_loss - loss_history(iter)) < tol && iter>1
        Y_abs = Y_abs_new; break;
    end
    Y_abs = Y_abs_new; prev_loss = loss_history(iter);
end

Y_abs_final = Y_abs;
Y_prop_final = Y_abs_final ./ sum(Y_abs_final,2);

%% —— Evaluation (Log Space) ——
X_fit_log = b_iter' + Y_abs_final*A_iter';
Res_raw   = X_log - X_fit_log;
Res_w     = Res_raw;
RMSE      = sqrt(mean(Res_w.^2,2,'omitnan'));

%% —— Similarity and Metrics —— 
fragNames = {'Short','Medium','Long'};
refP = [0.9 0.05 0.05];
brayCurtisSim = 1 - sum(abs(Y_prop_final - refP),2) ./ (sum(Y_prop_final,2) + sum(refP) + eps);
cosineSim = brayCurtisSim;

pearsonR = zeros(size(Y_prop_final,1),1);
for i = 1:size(Y_prop_final,1)
    a = log1p(refP); b = log1p(Y_prop_final(i,:));
    R = corrcoef(a(:), b(:)); pearsonR(i) = R(1,2);
end

%% —— Export CSVs (Consistent with KRAS Naming) ——
T_abs = array2table(Y_abs_final, 'VariableNames', fragNames); 
T_abs.Sample = sampleNames; T_abs.Type = normType; T_abs.GIT = normGIT;
writetable(T_abs, fullfile(outdir, 'final_absolute_fragment_concentration.csv'));

T_prop = array2table(Y_prop_final, 'VariableNames', fragNames); 
T_prop.Sample = sampleNames;
writetable(T_prop, fullfile(outdir, 'final_predicted_fragment_proportion.csv'));

T_fitAbs = array2table(Y_abs_final, 'VariableNames', fragNames); 
T_fitAbs.Sample = sampleNames;
writetable(T_fitAbs, fullfile(outdir, 'predicted_absolute_concentration.csv'));

T_r2 = table((1:length(R2_history))', R2_history, 'VariableNames', {'Iteration','R_squared'});
writetable(T_r2, fullfile(outdir, 'R2_history.csv'));

T_err = table(sampleNames, RMSE, brayCurtisSim, pearsonR, ...
    'VariableNames', {'Sample','RMSE_weighted_log','Similarity_BrayCurtis','Pearson_log1p'});
writetable(T_err, fullfile(outdir, 'error_metrics_per_sample.csv'));

%% —— Unified Plotting Style —— 
set(groot, 'defaultAxesFontSize', 30);
set(groot, 'defaultAxesLineWidth', 1.5);
set(groot, 'defaultLineLineWidth', 2);
set(groot, 'defaultBarLineWidth', 1.5);
set(groot, 'defaultAxesFontWeight', 'bold');
set(groot, 'defaultLineMarkerSize', 10);
colorMap = lines(3);
cmap = lines(numel(sampleNames));

%% Image Output (Consistent with KRAS)
% (1) Static R² Bar Chart
figure; bar(R2_history, 0.4, 'FaceColor', colorMap(1,:));
ylabel('R^2', 'FontWeight', 'bold'); xlabel('Iteration', 'FontWeight', 'bold'); title('R^2 per Iteration');
grid on; saveas(gcf, fullfile(outdir, 'iterative_R2_convergence.png'));

% (2) Predicted Absolute Concentration Grouped Bar Chart
figure;
b = bar(Y_abs_final, 1, 'grouped'); for k = 1:3, b(k).FaceColor = colorMap(k,:); end
title('Predicted Absolute Fragment Concentrations');
ylabel('EGFR EFIRM (cps/uL)', 'FontWeight', 'bold'); xlabel('Sample', 'FontWeight', 'bold');
set(gca, 'XTick', 1:numel(sampleNames), 'XTickLabel', sampleNames, 'FontWeight', 'bold'); xtickangle(45);
grid on; legend(fragNames, 'Location', 'northeast');
saveas(gcf, fullfile(outdir, 'fragment_absolute_barplot.png'));

% (3) Predicted Fragment Proportion Grouped Bar Chart
figure;
b = bar(Y_prop_final, 1, 'grouped'); for k = 1:3, b(k).FaceColor = colorMap(k,:); end
title('Predicted Fragment Proportions after Iteration');
ylabel('Proportion', 'FontWeight', 'bold'); xlabel('Sample', 'FontWeight', 'bold');
set(gca, 'XTick', 1:numel(sampleNames), 'XTickLabel', sampleNames, 'FontWeight', 'bold'); xtickangle(45);
grid on; legend(fragNames, 'Location', 'northeast');
saveas(gcf, fullfile(outdir, 'fragment_barplot_proportion.png'));

%% —— ADDED: Absolute Concentrations Grouped Bar Chart (usctDNA & mnctDNA only) ——
% usctDNA = Short = Y_abs_final(:,1)
% mnctDNA = Long  = Y_abs_final(:,3)
Y_us_mn = [Y_abs_final(:,1), Y_abs_final(:,3)];   % [usctDNA, mnctDNA]
N = size(Y_us_mn,1);

% Fixed Colors
col_us = [0.2 0.45 0.9];   % Blue
col_mn = [0.85 0.1 0.55];  % Magenta-pink

figure('Units','pixels','Position',[100,100,1600,600]);
b2 = bar(Y_us_mn, 1, 'grouped');
b2(1).FaceColor = col_us; 
b2(2).FaceColor = col_mn;
title('Predicted Absolute Fragment Concentrations (usctDNA & mnctDNA)');
ylabel('EGFR EFIRM (cps/\muL)', 'FontWeight', 'bold');
xlabel('NSCLC samples', 'FontWeight', 'bold');

% Keep axis lines and labels but hide X-axis tick numbers
set(gca, 'XTick', 1:N, 'XTickLabel', [], ...
    'XColor','k','FontWeight','bold');
xlabel('NSCLC samples','FontWeight','bold'); 
xtickangle(0); grid on;
legend({'usctDNA','mnctDNA'}, 'Location', 'northeast');
saveas(gcf, fullfile(outdir, 'fragment_absolute_barplot_usctDNA_mnctDNA.png'));

%% —— ADDED: Box + Scatter Plots for usctDNA / mnctDNA Classes ——
% Description:
% - usctDNA = Short; mnctDNA = Long
% - Uses only column 1 (Short) and 3 (Long) of Y_abs_final
% - Exports a CSV with only these two classes for reporting

us_mn_names = {'usctDNA','mnctDNA'};
Y_us_mn     = [Y_abs_final(:,1), Y_abs_final(:,3)];  % [Short, Long] → [usctDNA, mnctDNA]

% 1) Plot boxchart + scatter
Y_us_mn_vec   = reshape(Y_us_mn', [], 1);                            % Vectorization
Group_us_mn   = repmat(us_mn_names', numel(sampleNames), 1);         % Group labels
fragCats_usmn = categorical(Group_us_mn, us_mn_names, 'Ordinal',true);

figure('Units','pixels','Position',[120,120,800,700]); hold on;
boxchart(fragCats_usmn, Y_us_mn_vec, 'BoxWidth',0.45, 'BoxFaceColor',[0.65 0.65 0.90]);

% Overlay sample points with jitter (colored by sample)
cmap_usmn = lines(numel(sampleNames));
for i = 1:numel(sampleNames)
    % usctDNA (Original: Short)
    x1 = 1 + 0.12*(rand-0.5);
    scatter(x1, Y_abs_final(i,1), 60, 'filled', 'MarkerFaceColor', cmap_usmn(i,:));
    % mnctDNA (Original: Long)
    x2 = 2 + 0.12*(rand-0.5);
    scatter(x2, Y_abs_final(i,3), 60, 'filled', 'MarkerFaceColor', cmap_usmn(i,:));
end

title('Absolute Concentrations: usctDNA & mnctDNA');
ylabel('EGFR EFIRM (cps/\muL)','FontWeight','bold'); 
xlabel('Fragment Class','FontWeight','bold');
grid on; 
saveas(gcf, fullfile(outdir, 'box_usctDNA_mnctDNA_absolute.png'));

% 2) Export CSV with usctDNA / mnctDNA only (cps/uL)
T_us_mn = table(sampleNames, normType, normGIT, ...
    Y_abs_final(:,1), Y_abs_final(:,3), ...
    'VariableNames', {'Sample','Type','GIT','usctDNA','mnctDNA'});
writetable(T_us_mn, fullfile(outdir, 'absolute_usctDNA_mnctDNA_only.csv'));

%% —— Plasma / Saliva Dual Plots (usctDNA & mnctDNA box+scatter) ——
typeOrder   = {'P','S'};
typeTitles  = containers.Map({'P','S'}, {'Plasma (P)','Saliva (S)'});

% Color definitions
col_us = [0.2 0.45 0.9];   % usctDNA Blue
col_mn = [0.85 0.1 0.55];  % mnctDNA Magenta-pink

% Determine global Y-axis range for consistency across Plasma and Saliva
ymax = max([Y_abs_final(normType=="P",1); Y_abs_final(normType=="P",3); ...
            Y_abs_final(normType=="S",1); Y_abs_final(normType=="S",3)], [], 'omitnan');
ymax = ceil(ymax*1.1);  % Add safety margin

fig = figure('Units','pixels','Position',[120,120,1600,800]); % 2:1 aspect ratio
tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');

for ti = 1:numel(typeOrder)
    Tcode = typeOrder{ti};
    idx   = normType == Tcode;
    if ~any(idx)
        nexttile; axis off; title([typeTitles(Tcode) ' (no data)']); continue;
    end
    
    % Get absolute concentrations
    Y_us = Y_abs_final(idx,1);
    Y_mn = Y_abs_final(idx,3);
    
    % Assemble boxchart inputs
    Y_vec = [Y_us; Y_mn];
    G_vec = [repmat(categorical({'usctDNA'}), numel(Y_us),1); ...
             repmat(categorical({'mnctDNA'}), numel(Y_mn),1)];
    
    ax = nexttile; hold(ax,'on');
    
    % Boxcharts with distinct colors
    sel_us = (G_vec == 'usctDNA');
    sel_mn = (G_vec == 'mnctDNA');
    boxchart(G_vec(sel_us), Y_vec(sel_us), 'BoxWidth',0.5, 'BoxFaceColor', col_us, 'MarkerStyle','none');
    boxchart(G_vec(sel_mn), Y_vec(sel_mn), 'BoxWidth',0.5, 'BoxFaceColor', col_mn, 'MarkerStyle','none');
    
    % Overlay scatter points with jitter
    rng(2);
    x1 = 1 + 0.12*(rand(numel(Y_us),1)-0.5);
    x2 = 2 + 0.12*(rand(numel(Y_mn),1)-0.5);
    scatter(x1, Y_us, 30, 'filled', 'MarkerFaceColor', col_us, 'MarkerEdgeColor','k','LineWidth',0.3);
    scatter(x2, Y_mn, 30, 'filled', 'MarkerFaceColor', col_mn, 'MarkerEdgeColor','k','LineWidth',0.3);
    
    title(typeTitles(Tcode));
    ylabel('EGFR EFIRM (cps/\muL)','FontWeight','bold');
    ylim([0 ymax]);  % Unified Y-axis range
    grid on;
    set(ax,'FontWeight','bold');
end

saveas(fig, fullfile(outdir, 'box_usctDNA_mnctDNA_Plasma_Saliva.png'));

%% —— Evaluation Charts (Bar and Histogram) ——
% Cosine Similarity bar chart (Uses Bray-Curtis but compatible with legacy naming)
figure; bar(cosineSim, 0.4, 'FaceColor', colorMap(3,:));
title('Cosine Similarity per Sample (Absolute)'); ylabel('Cosine Similarity', 'FontWeight', 'bold'); xlabel('Sample', 'FontWeight', 'bold');
set(gca, 'XTick', 1:numel(sampleNames), 'XTickLabel', sampleNames, 'FontWeight', 'bold'); xtickangle(45);
grid on; saveas(gcf, fullfile(outdir, 'bar_CosineSimilarity_absolute.png'));

% RMSE Histogram and QQ Plot
edges = 0:0.02:max(max(RMSE)+0.1, 0.5);
figure; histogram(RMSE, edges, 'FaceColor', colorMap(2,:));
title('Histogram of RMSE (Absolute)'); xlabel('RMSE', 'FontWeight', 'bold'); ylabel('Frequency', 'FontWeight', 'bold');
grid on; saveas(gcf, fullfile(outdir, 'histogram_RMSE_absolute.png'));

figure('Units','pixels','Position',[100,100,900,900]); qqplot(RMSE(:));
title('QQ Plot of RMSE (Absolute)'); grid on; saveas(gcf, fullfile(outdir, 'qqplot_RMSE_absolute.png'));

% Cosine Similarity Histogram and QQ Plot
edges2 = 0:0.02:1;
figure; histogram(cosineSim, edges2, 'FaceColor', colorMap(3,:));
title('Histogram of Cosine Similarity (Absolute)'); xlabel('Cosine Similarity','FontWeight','bold'); ylabel('Frequency','FontWeight','bold');
grid on; saveas(gcf, fullfile(outdir, 'histogram_CosineSimilarity_absolute.png'));

figure('Units','pixels','Position',[100,100,900,900]); qqplot(cosineSim(:));
title('QQ Plot of Cosine Similarity (Absolute)'); grid on;
saveas(gcf, fullfile(outdir, 'qqplot_CosineSimilarity_absolute.png'));

% Bray-Curtis QQ Plot
figure('Units','pixels','Position',[100,100,900,900]); qqplot(brayCurtisSim(:));
title('QQ Plot of Bray-Curtis Similarity (Absolute)'); grid on;
saveas(gcf, fullfile(outdir, 'qqplot_BrayCurtis_absolute.png'));

% Residuals Histogram with Normal and KDE fitting
fig = figure('Units','pixels','Position',[100,100,900,600]); hold on;
h = histogram(Res_w(:),30,'Normalization','pdf','FaceColor',[0.85 0.4 0.2],'EdgeColor','none');
mu = mean(Res_w(:)); sg = std(Res_w(:));
xx = linspace(min(Res_w(:)), max(Res_w(:)), 200);
plot(xx, normpdf(xx,mu,sg),'b-','LineWidth',2);
[f,xi] = ksdensity(Res_w(:)); plot(xi,f,'k--','LineWidth',2);
title('Histogram of Weighted Residuals','FontSize',22,'FontWeight','bold');
xlabel('Weighted residuals in log1p space','FontWeight','bold'); ylabel('Density','FontWeight','bold');
legend({'Histogram (pdf)','Normal fit','KDE fit'},'Location','best'); grid on;
saveas(fig, fullfile(outdir,'hist_weighted_residuals_with_fit.png'));

% Scatter: Residuals vs Predicted
fig = figure('Units','pixels','Position',[100,100,900,600]);
scatter(X_fit_log(:), Res_w(:), 36, 'filled'); grid on; yline(0,'--');
xlabel('Predicted (log1p space)','FontWeight','bold'); ylabel('Weighted residual','FontWeight','bold');
title('Residuals vs Predicted','FontSize',22,'FontWeight','bold');
saveas(fig, fullfile(outdir,'scatter_residual_vs_pred.png'));

%% —— Stacked Bar Chart + ddPCR Overlay (if available) —— 
if any(~isnan(ddPCR_vec))
    fig = figure('Units','pixels','Position',[100,100,1100,700]); hold on; box on; grid on;
    barH = bar([Y_abs_final(:,1) Y_abs_final(:,2) Y_abs_final(:,3)], 'stacked', 'BarWidth',0.75);
    barH(1).FaceColor=[0.30 0.50 0.85]; barH(2).FaceColor=[0.95 0.55 0.35]; barH(3).FaceColor=[0.35 0.75 0.55];
    for k=1:3, barH(k).EdgeColor='none'; end
    plot(1:numel(ddPCR_vec), ddPCR_vec, '-ok','LineWidth',2,'MarkerFaceColor','k');
    set(gca,'XTick',1:numel(sampleNames),'XTickLabel',sampleNames,'XTickLabelRotation',30,'FontWeight','bold');
    ylabel('Concentration (cps/\muL)','FontWeight','bold');
    title('Stacked Fragment Concentration with ddPCR','FontSize',22,'FontWeight','bold');
    legend({'Short','Medium','Long','ddPCR'},'Location','northwest'); 
    saveas(fig, fullfile(outdir,'stacked_fragments_with_ddPCR.png'));
end

%% === Parallel-coordinates style plot (Short–Medium–Long) ===
fig = figure('Units','pixels','Position',[120,120,900,600]); hold on;
axis([0.8 3.2 0 1]); box on; grid on; set(gca,'Layer','top');
set(gca,'XTick',[1 2 3],'XTickLabel',{'Short %','Medium %','Long %'}, ...
    'FontWeight','bold','FontSize',16);

% Plot three vertical "axes" lines
plot([1 1],[0 1],'k-','LineWidth',1);
plot([2 2],[0 1],'k-','LineWidth',1);
plot([3 3],[0 1],'k-','LineWidth',1);

% Color based on Short proportion
shortVals = Y_prop_final(:,1);
cmap = parula(256);
cidx = max(1, min(256, round( 1 + 255 * (shortVals - min(shortVals)) / max(eps, (max(shortVals)-min(shortVals))) )));

% Draw trajectories using smoothed Bezier curves
for i = 1:size(Y_prop_final,1)
    s = Y_prop_final(i,1);
    m = Y_prop_final(i,2);
    l = Y_prop_final(i,3);
    col = cmap(cidx(i),:);
    
    % Section 1: (1,s) -> (2,m)
    x1 = 1; y1 = s;
    x2 = 2; y2 = m;
    cx = (x1+x2)/2;
    cy = (y1+y2)/2 + 0.25*(y1 - y2); % Curvature factor
    t  = linspace(0,1,60);
    bx1 = (1-t).^2*x1 + 2*(1-t).*t*cx + t.^2*x2;
    by1 = (1-t).^2*y1 + 2*(1-t).*t*cy + t.^2*y2;
    
    % Section 2: (2,m) -> (3,l)
    x3 = 3; y3 = l;
    cx2 = (x2+x3)/2;
    cy2 = (y2+y3)/2 + 0.25*(y3 - y2);
    bx2 = (1-t).^2*x2 + 2*(1-t).*t*cx2 + t.^2*x3;
    by2 = (1-t).^2*y2 + 2*(1-t).*t*cy2 + t.^2*y3;
    
    plot(bx1, by1, '-', 'Color', col, 'LineWidth', 1.5);
    plot(bx2, by2, '-', 'Color', col, 'LineWidth', 1.5);
end

ylim([0 1]); yticks(0:0.1:1);
ylabel('Proportion', 'FontWeight','bold');

% Colorbar showing Short proportion
cb = colorbar; cb.Label.String = 'Short %';
cb.Label.FontWeight = 'bold';
if max(shortVals) > min(shortVals)
    ticks = linspace(min(shortVals), max(shortVals), 7);
else
    ticks = shortVals(1);
end
cb.Ticks = (ticks - min(shortVals)) / max(eps, (max(shortVals)-min(shortVals)));
cb.TickLabels = compose('%.3f', ticks);

title('Fragment Proportion Trajectories (Short–Medium–Long)', ...
      'FontSize',18,'FontWeight','bold');
saveas(fig, fullfile(outdir, 'parallel_fragments_short_medium_long.png'));

fprintf('Analysis complete. Results saved in: %s\n', outdir);
clearvars -except outdir

%% —— Internal Functions —— 
function namesOut = make_unique_names(namesIn)
    % Resolves duplicate names by adding _2, _3 suffixes
    namesOut = strings(size(namesIn));
    mp = containers.Map('KeyType','char','ValueType','int32');
    for i=1:numel(namesIn)
        key = char(namesIn(i));
        if isKey(mp,key)
            mp(key) = mp(key)+1;
            namesOut(i) = string(key) + "_" + string(mp(key));
        else
            mp(key) = 1;
            namesOut(i) = string(key);
        end
    end
end

function f = obj_huber_p(p, Aiter, biter, xlog, w, ktotal, delta, gamma)
    % Huber loss objective function with entropy regularization for proportions p
    pred = biter + Aiter*(ktotal * p(:));
    r = xlog - pred;
    idx = abs(r) <= delta;
    hub = zeros(size(r));
    hub(idx)  = 0.5 * r(idx).^2;
    hub(~idx) = delta * (abs(r(~idx)) - 0.5*delta);
    
    % Entropy penalty: sum(p log p) to avoid edge solutions
    pen = gamma * sum(p(:) .* log(max(p(:), 1e-12)));
    f = sum(w(:).*hub(:)) + pen;
end