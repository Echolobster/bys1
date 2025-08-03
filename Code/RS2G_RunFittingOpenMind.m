%% RS2G_fittingOpenMind
%
%   Fits the BLSbiased model to a set of subjects and (optional) saves
%   results.
%
%%

% Gobal Variables
trialWin = [100 Inf];
interval_N = 1:2;
MinMax = [0 1000];
outlier = 5;
FitAll = 1;
FitRuns = 0;
SaveFlg = 1;

% Fit parameters
fparams.fittype = 'BLSbiasedLapse';
fparams.method = 'quad';
fparams.init = 'estb';
fparams.dx = 10;
fparams.trialtypes = 1;%[1 2];

% Bias/Variance bootstrap parameters
bootparams.nbootstraps = 100;
bootparams.nsamps = 50;

% Parameters for calculating expected aim times
TAexpectation.method = 'numerical';
TAexpectation.trialtypes = [1 2];
TAexpectation.ts_vec = (550:10:1050)';
TAexpectation.simtrials = 10000;

% Subject list
Subjects = {'LB'};%{'CM','CV','GB','LB','PG','SC','TA','VD','VR'};

% Run through each subject and fit
for SubjectN = 1:length(Subjects)
    disp(['Subject ' Subjects{SubjectN}])
    
    % Load the data
    d = load([Subjects{SubjectN} '_RS2G_psychophysics']);
    runs = 2:d.runs;
    
    
    if SaveFlg
        SaveFileBase = ['/home/swegger/Projects/RS2G_psychophysics/' Subjects{SubjectN} '/' Subjects{SubjectN} '_BLSbiasedFitResults'];
        SaveParam = [SaveFileBase datestr(now,'yyyymmdd')];
    else
        SaveParam = 'No';
    end
    
    
    % Fit all the data
    if FitAll
        [mtp, stdtp, bias, variance, rmse, wm, wp, b, pval, weber, tsIn, tpIn, trialsIn, ts_in, tp_in, Trials_sorted, estb, ta, G, Llikelihood, BIASs, VARs, lapse, lapseTrials] = RS2G_psychophysicsPoolAnalysis(d,'runs',runs(2:end),'interval_N',interval_N,'Fit',fparams,'Plot','No',...
            'outlier',outlier,'ConflictType','equal','MinMax',MinMax,'Bootstrap',bootparams,...
            'TAexpectation',TAexpectation,'trialWin',trialWin,'Save',SaveParam,...
            'interval_N',interval_N);
    end
    
    % Fit each run
    if FitRuns
        for i = runs
            [mtp_r{i}, stdtp_r{i}, bias_r{i}, variance_r{i}, rmse_r{i}, wm_r{i}, wp_r{i}, b_r{i}, pval_r{i}, weber_r{i}, tsIn_r{i}, tpIn_r{i}, trialsIn_r{i}, ts_in_r{i}, tp_in_r{i}, Trials_sorted_r{i}, estb_r{i}, ta_r{i}, G_r{i}, Llikelihood_r{i}, BIASs_r{i}, VARs_r{i}, lapse_r{i}, lapseTrials_r{i}] = RS2G_psychophysicsPoolAnalysis(d,...
                'runs',runs,'interval_N',interval_N,'Fit',fparams,'Plot','No',...
                'outlier',outlier,'ConflictType','equal','MinMax',MinMax,'Bootstrap',bootparams,...
                'TAexpectation',TAexpectation,'trialWin',trialWin);
            
        end
    end
        
    
end
