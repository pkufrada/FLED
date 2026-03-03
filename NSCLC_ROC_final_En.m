%% -----------------------------------------------------------
%  EFIRM ROC Analysis Script
%  Automatically creates a timestamped folder to save all data and plots.
%  Updates:
%   1) Increased overall plot font sizes by +2.
%   2) Saves the current script as a .txt file for record-keeping.
%   3) Cleaned up redundant code segments and isolated characters.
% ------------------------------------------------------------

%% Input Filenames (Modify as needed)
plasmaFile = '20250829 Plasma GIT-NEG.xlsx';
salivaFile = '20250829 Saliva GIT-NEG.xlsx';

% Mapping: short = usctDNA, long = mnctDNA
shortCol = 1;       % Column A
longCol  = 3;       % Column C
labelCol = 7;       % Column G (tissue genotype)

%% -----------------------------------------------------------
% 0. Font Size Configuration (Increased by +2 from baseline)
%% -----------------------------------------------------------
FS_AX    = 16;   % Axis ticks (was 14)
FS_LAB   = 18;   % Axis labels (was 16)
FS_TITLE = 22;   % Plot titles (was 20)
FS_AUC   = 16;   % AUC text (was 14)

%% -----------------------------------------------------------
% 1. Create Output Directory with Current Timestamp
%% -----------------------------------------------------------
timestamp = datestr(now,'yyyymmdd_HHMMSS');
outdir = ['output_' timestamp];
mkdir(outdir);
fprintf('Output directory created: %s\n', outdir);

%% -----------------------------------------------------------
% 2. Execute ROC Computation
%% -----------------------------------------------------------
rocPlasma = computeROC(plasmaFile, shortCol, longCol, labelCol, outdir, "Plasma");
rocSaliva = computeROC(salivaFile, shortCol, longCol, labelCol, outdir, "Saliva");

%% Save Workspace Results
save(fullfile(outdir, 'ROC_results.mat'), 'rocPlasma', 'rocSaliva');

%% -----------------------------------------------------------
% 3. Export ROC Data Tables (Excel) - Padding with NaN for length consistency
%% -----------------------------------------------------------
% Helper: Pad vector v to length n using NaNs
padcol = @(v,n) [v(:); nan(n - numel(v), 1)];

% ----- Plasma Data Processing -----
nP = max([ numel(rocPlasma.fpr_short), ...
           numel(rocPlasma.tpr_short), ...
           numel(rocPlasma.fpr_long), ...
           numel(rocPlasma.tpr_long) ]);
P_fpr_short = padcol(rocPlasma.fpr_short, nP);
P_tpr_short = padcol(rocPlasma.tpr_short, nP);
P_fpr_long  = padcol(rocPlasma.fpr_long,  nP);
P_tpr_long  = padcol(rocPlasma.tpr_long,  nP);

ROCtablePlasma = table(P_fpr_short, P_tpr_short, P_fpr_long, P_tpr_long, ...
    'VariableNames', {'Plasma_FPR_short','Plasma_TPR_short', ...
                      'Plasma_FPR_long','Plasma_TPR_long'});

% ----- Saliva Data Processing -----
nS = max([ numel(rocSaliva.fpr_short), ...
           numel(rocSaliva.tpr_short), ...
           numel(rocSaliva.fpr_long), ...
           numel(rocSaliva.tpr_long) ]);
S_fpr_short = padcol(rocSaliva.fpr_short, nS);
S_tpr_short = padcol(rocSaliva.tpr_short, nS);
S_fpr_long  = padcol(rocSaliva.fpr_long,  nS);
S_tpr_long  = padcol(rocSaliva.tpr_long,  nS);

ROCtableSaliva = table(S_fpr_short, S_tpr_short, S_fpr_long, S_tpr_long, ...
    'VariableNames', {'Saliva_FPR_short','Saliva_TPR_short', ...
                      'Saliva_FPR_long','Saliva_TPR_long'});

% Export to a single Excel file with multiple sheets
outExcel = fullfile(outdir, 'ROC_data.xlsx');
writetable(ROCtablePlasma, outExcel, 'Sheet','Plasma');
writetable(ROCtableSaliva, outExcel, 'Sheet','Saliva');

%% -----------------------------------------------------------
% 4. Visualization: Combined ROC Plots (Plasma & Saliva)
%% -----------------------------------------------------------
fig = figure('Units','pixels', ...
             'Position',[200 100 900 1200], ...   
             'Name','NSCLC ROC Analysis', ...
             'NumberTitle','off');

% Color Definitions
usctColor = [0.56 0.76 0.96];   % usctDNA (short)
mnctColor = [0.96 0.68 0.84];   % mnctDNA (long)
diagColor = [0.6 0.6 0.6];      % Diagonal reference line

% Layout: 1x2 grid with tight spacing
tiledlayout(1, 2, 'TileSpacing','compact', 'Padding','tight');

%% ================= Subplot: Plasma =================
ax1 = nexttile;
hold(ax1, 'on'); box(ax1, 'on');

plot(ax1, rocPlasma.fpr_short, rocPlasma.tpr_short, '-', ...
    'LineWidth', 3, 'Color', usctColor);
plot(ax1, rocPlasma.fpr_long,  rocPlasma.tpr_long,  '-', ...
    'LineWidth', 3, 'Color', mnctColor);
plot(ax1, [0 1], [0 1], '--', 'Color', diagColor, 'LineWidth', 1.5);

grid(ax1, 'on');
ax1.GridLineStyle = '--';
ax1.GridColor     = [0.7 0.7 0.7];
ax1.GridAlpha     = 0.4;

xlim(ax1, [0 1]);
ylim(ax1, [0 1]);
axis(ax1, 'square');

set(ax1, 'FontName','Arial Narrow', ...
         'FontSize', FS_AX, ...
         'LineWidth', 1.6);

xlabel(ax1, '1 - Specificity (False Positive Rate)', ...
       'FontName','Arial Narrow', 'FontSize', FS_LAB);
ylabel(ax1, 'Sensitivity (True Positive Rate)', ...
       'FontName','Arial Narrow', 'FontSize', FS_LAB);
title(ax1, 'Plasma', ...
      'FontName','Arial Narrow', 'FontSize', FS_TITLE, 'FontWeight','bold');

% Annotate AUC Values
text(ax1, 0.35, 0.20, ...
    sprintf('usctDNA  AUC = %.2f', rocPlasma.auc_short), ...
    'Color', usctColor, ...
    'FontSize', FS_AUC, ...
    'FontName','Arial Narrow', ...
    'FontWeight','bold', ...
    'HorizontalAlignment','left', ...
    'Clipping','on');
text(ax1, 0.35, 0.12, ...
    sprintf('mnctDNA  AUC = %.2f', rocPlasma.auc_long), ...
    'Color', mnctColor, ...
    'FontSize', FS_AUC, ...
    'FontName','Arial Narrow', ...
    'FontWeight','bold', ...
    'HorizontalAlignment','left', ...
    'Clipping','on');
hold(ax1, 'off');

%% ================= Subplot: Saliva =================
ax2 = nexttile;
hold(ax2, 'on'); box(ax2, 'on');

plot(ax2, rocSaliva.fpr_short, rocSaliva.tpr_short, '-', ...
    'LineWidth', 3, 'Color', usctColor);
plot(ax2, rocSaliva.fpr_long,  rocSaliva.tpr_long,  '-', ...
    'LineWidth', 3, 'Color', mnctColor);
plot(ax2, [0 1], [0 1], '--', 'Color', diagColor, 'LineWidth', 1.5);

grid(ax2, 'on');
ax2.GridLineStyle = '--';
ax2.GridColor     = [0.7 0.7 0.7];
ax2.GridAlpha     = 0.4;

xlim(ax2, [0 1]);
ylim(ax2, [0 1]);
axis(ax2, 'square');

set(ax2, 'FontName','Arial Narrow', ...
         'FontSize', FS_AX, ...
         'LineWidth', 1.6);

xlabel(ax2, '1 - Specificity (False Positive Rate)', ...
       'FontName','Arial Narrow', 'FontSize', FS_LAB);
ylabel(ax2, 'Sensitivity (True Positive Rate)', ...
       'FontName','Arial Narrow', 'FontSize', FS_LAB);
title(ax2, 'Saliva', ...
      'FontName','Arial Narrow', 'FontSize', FS_TITLE, 'FontWeight','bold');

% Annotate AUC Values
text(ax2, 0.35, 0.20, ...
    sprintf('usctDNA  AUC = %.2f', rocSaliva.auc_short), ...
    'Color', usctColor, ...
    'FontSize', FS_AUC, ...
    'FontName','Arial Narrow', ...
    'FontWeight','bold', ...
    'HorizontalAlignment','left', ...
    'Clipping','on');
text(ax2, 0.35, 0.12, ...
    sprintf('mnctDNA  AUC = %.2f', rocSaliva.auc_long), ...
    'Color', mnctColor, ...
    'FontSize', FS_AUC, ...
    'FontName','Arial Narrow', ...
    'FontWeight','bold', ...
    'HorizontalAlignment','left', ...
    'Clipping','on');
hold(ax2, 'off');

%% -----------------------------------------------------------
% 5. Save Figures
%% -----------------------------------------------------------
saveas(fig, fullfile(outdir, 'ROC_combined.png'));
saveas(fig, fullfile(outdir, 'ROC_combined.pdf'));

%% -----------------------------------------------------------
% 6. Backup Current Script for Reproducibility (.m and .txt)
%% -----------------------------------------------------------
stackInfo = dbstack('-completenames');
if ~isempty(stackInfo)
    currentScriptFull = stackInfo(1).file;   
else
    currentScriptFull = mfilename('fullpath');
end
[scriptPath, scriptName, ~] = fileparts(currentScriptFull);

% 6.1 Save as .m file
try
    copyfile(fullfile(scriptPath, [scriptName '.m']), ...
             fullfile(outdir, [scriptName '_saved.m']));
    fprintf('Code backed up as: %s\n', fullfile(outdir, [scriptName '_saved.m']));
catch ME
    warning('Could not back up .m script: %s', ME.message);
end

% 6.2 Save as .txt file for universal reading
try
    codeStr = fileread(currentScriptFull);
    txtPath = fullfile(outdir, [scriptName '_saved.txt']);
    fid = fopen(txtPath, 'w');
    fwrite(fid, codeStr);
    fclose(fid);
    fprintf('Code backed up as: %s\n', txtPath);
catch ME
    warning('Could not back up as .txt: %s', ME.message);
end

fprintf('All data and charts saved to directory: %s\n', outdir);

%% ===========================================================
%  Internal Function: computeROC
%% ===========================================================
function rocData = computeROC(fname, shortCol, longCol, labelCol, outdir, prefix)
    T = readtable(fname,'FileType','spreadsheet');
    
    shortScore = T{:, shortCol};      % usctDNA scores
    longScore  = T{:, longCol};       % mnctDNA scores
    gtRaw      = string(T{:, labelCol});
    labels     = gtRaw == "L858R";    % Define positive class
    
    % Save data extraction check (for debugging)
    debugTable = table(shortScore, longScore, labels, gtRaw, ...
        'VariableNames', {'Short_usct','Long_mnct','Label_binary','Label_raw'});
    writetable(debugTable, fullfile(outdir, prefix + "_RawDataCheck.xlsx"));
    
    % ROC Calculations
    [fpr_s, tpr_s, ~, auc_s] = perfcurve(labels, shortScore, 1, ...
        'XCrit','FPR','YCrit','TPR');
    [fpr_l, tpr_l, ~, auc_l] = perfcurve(labels, longScore, 1, ...
        'XCrit','FPR','YCrit','TPR');
    
    rocData.fpr_short = fpr_s;
    rocData.tpr_short = tpr_s;
    rocData.auc_short = auc_s;
    
    rocData.fpr_long  = fpr_l;
    rocData.tpr_long  = tpr_l;
    rocData.auc_long  = auc_l;
    
    % Store raw scores and binary labels
    rocData.shortScore = shortScore;
    rocData.longScore  = longScore;
    rocData.labels     = labels;
end