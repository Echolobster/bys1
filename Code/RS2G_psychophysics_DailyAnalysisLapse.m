function d = RS2G_psychophysics_DailyAnalysisLapse(sname)
%% RS2G_psychophysics Daily Analysis
%
%   Extracts and fits a subjects data; saves fit results to data structure
%
%%

% Variables
trialWin = [100 Inf];
interval_N = 1:2;
MinMax = [0 1000];
outlier = Inf;
SaveLocation = 'LocalHeader';

% Fit parameters
fparams.fittype = 'BLSbiasedLapse';
fparams.method = 'quad';
fparams.init = 'estb';
fparams.dx = 10;
fparams.trialtypes = [1 2];

% Bias/Variance bootstrap parameters
bootparams.nbootstraps = 100;
bootparams.nsamps = 50;

% Parameters for calculating expected aim times
TAexpectation.method = 'numerical';
TAexpectation.trialtypes = [1 2];
TAexpectation.ts_vec = (550:10:1050)';
TAexpectation.simtrials = 10000;

% Load data structure and extract the data
d = SetupMWorksData('RS2G_psychophysics',sname,'SaveLocation','LocalHeader');


% Determine which runs have not been fit
fitflg = nan(d.runs,1);
if ~isfield(d,'modelParams')
    d = addFieldToProj('RS2G_psychophysics',sname,'modelParams','d',d,'SaveLocation','LocalHeader');
end
for i = 1:d.runs
    if length(d.modelParams) < i || ~isfield(d.modelParams{i},'lapse')
        d.modelParams{i} = [];
        fitflg(i) = 1;
    else
        fitflg(i) = isempty(d.modelParams{i});
    end
end

% Fit the data and plot results
for i = 1:d.runs
    if fitflg(i)
        [mtp, stdtp, bias, variance, rmse, wm, wp, b, pval, weber, tsIn, tpIn, trialsIn, ts_in, tp_in, Trials_sorted, estb, ta, G, Llikelihood, BIASs, VARs, lapse, lapseTrials, simbias, simv] = RS2G_psychophysicsPoolAnalysis(d,'runs',i,'interval_N',interval_N,'Fit',fparams,'Plot','All',...
                                                                                                                                                                                    'outlier',outlier,'ConflictType','equal','MinMax',MinMax,'Bootstrap',bootparams,...
                                                                                                                                                                                    'TAexpectation',TAexpectation,'trialWin',trialWin);
        d.modelParams{i}.wm = wm;
        d.modelParams{i}.wp = wp;
        d.modelParams{i}.b = b;
        d.modelParams{i}.lapse = lapse;
        
        % Save the results
        switch SaveLocation
            case 'default'
                save(d.datafilepath, '-struct', 'd','-v7.3');
            case 'LocalHeader'
                A = fileread('/usr/local/matlab/HeaderFile');
                word = regexp(A,'savepath = .*;','match');
                savepath = [word{1}(13:end-2) '/' d.projname];
                filename = sprintf('%s_%s.mat', d.sname, d.projname);
                save([savepath '/' d.sname '/' filename], '-struct', 'd', '-v7.3');
        end
        
    end
end


