function RS2G_fittingOpenMind(SubjectN,varargin)
%% RS2G_fittingOpenMind
%
%   Fits the BLSbiased model to a set of subjects and (optional) saves
%   results.
%
%%

%% Defaults

%% Gobal Variables
Subjects = {'CM','CV','GB','LB','PG','RM','SC','SM','SWE','TA','VD','VR'};

trialWin = [100 Inf];
interval_N = 1:2;
MinMax = [0 1000];
outlier = Inf;
FitAll = 1;
FitRuns = 1;

% Fit parameters
fparams_default.fittype = {'SubOptMemBias','BLSbiasedLapse'};%,'WeightedLinear'};     % Specifies which models to fit to the data
fparams_default.modelUsed = 1;                                      % Specifies which model to use in subsequent model-based analysis
fparams_default.method = 'quad';                                    % Integration method in model fitting
fparams_default.init = 'estb';                                      % Initialization values/method for model fitting
fparams_default.dx = 10;                                            % Step size of integration (in ms)
fparams_default.trialtypes = [1 2];                                 % Trial types to fit
fparams_default.ObsAct = 0;

% Cross validation control
fparams_default.CrossValidation.Type = 'LNOCV';                              % Cross validation method
fparams_default.CrossValidation.N = 100;                                       % Left out trials for validation

% Model evidence control
fparams_default.ModelEvidence.method = 'none';
fparams_defaults.ModelEvidence.paramLimits = [0.0001 1;...
                             0.0001 1;...
                             -200 200;...
                             0 1];
fparams_default.ModelEvidence.integrationMethod = 'quad';
fparams_default.ModelEvidence.integrationOptions.dx = 0.1;
fparams_default.ModelEvidence.OpenMind = 1;

% Bias/Variance bootstrap parameters
bootparams.nbootstraps = 100;
bootparams.nsamps = 500;

% Parameters for calculating expected aim times
TAexpectation.method = 'numerical';
TAexpectation.trialtypes = [1 2];
TAexpectation.ts_vec = (550:10:1050)';
TAexpectation.simtrials = 10000;

SaveOpts_default.SaveFlg = true;

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'SubjectN')
addParameter(Parser,'SaveOpts',SaveOpts_default)
addParameter(Parser,'fparams',fparams_default)
addParameter(Parser,'TestSave',false)

parse(Parser,SubjectN,varargin{:})

SubjectN = Parser.Results.SubjectN;
SaveOpts = Parser.Results.SaveOpts;
fparams = Parser.Results.fparams;
TestSave = Parser.Results.TestSave;

fparams = asign_fparams(fparams);
if SaveOpts.SaveFlg && ~isfield(SaveOpts,'SaveFile')
    SaveOpts.SaveFile = ['/om/user/swegger/Projects/RS2G_psychophysics/'...
        Subjects{SubjectN} '/' Subjects{SubjectN} '_' ...
        fparams.fittype{fparams.modelUsed} '_ObsAct' num2str(fparams.ObsAct)...
        '_' datestr(now,'yyyymmdd')];
end
%% Run through each subject and fit
disp(['Subject ' Subjects{SubjectN}])

% Load the data
d = load([Subjects{SubjectN} '_RS2G_psychophysics']);
runs = 2:d.runs;


if SaveOpts.SaveFlg
    SaveParam = SaveOpts.SaveFile;
else
    SaveParam = 'No';
end


% Fit all the data
if FitAll && ~TestSave
    [mtp, stdtp, bias, variance, rmse, wm, wp, b, pval, weber, tsIn, tpIn,...
        trialsIn, ts_in, tp_in, Trials_sorted, estb, ta, G, Llikelihood,...
        BIASs, VARs, lapse, lapseTrials] = RS2G_psychophysicsPoolAnalysis(...
        d,'runs',runs,'interval_N',interval_N,'Fit',fparams,'Plot','No',...
        'outlier',outlier,'ConflictType','equal','MinMax',MinMax,'Bootstrap',bootparams,...
        'TAexpectation',TAexpectation,'trialWin',trialWin,'Save',SaveParam);
    
elseif TestSave && SaveOpts.SaveFlg
    save(SaveOpts.SaveFile)
end

%% Functions

function fparams = asign_fparams(fparams)
if ~isfield(fparams,'fittype')
    fparams.fittype = {'SubOptMemBias','BLSbiasedLapse'};
end
if ~isfield(fparams,'modelUsed')
    fparams.modelUsed = 1;                                      % Specifies which model to use in subsequent model-based analysis
end
if ~isfield(fparams,'method')
    fparams.method = 'quad';                                    % Integration method in model fitting
end
if ~isfield(fparams,'init')
    fparams.init = 'estb';                                      % Initialization values/method for model fitting
end
if ~isfield(fparams,'dx')
    fparams.dx = 10;                                            % Step size of integration (in ms)
end
if ~isfield(fparams,'trialtypes')
    fparams.trialtypes = [1 2];                                 % Trial types to fit
end
if ~isfield(fparams,'ObsAct')
    fparams.ObsAct = 0;
end

% Cross validation control
if ~isfield(fparams,'CrossValidation')
    fparams.CrossValidation.Type = 'LNOCV';                              % Cross validation method
    fparams.CrossValidation.N = 100;                                       % Left out trials for validation
end

% Model evidence control
if ~isfield(fparams,'ModelEvidence')
    fparams.ModelEvidence.method = 'none';
    fparams.ModelEvidence.paramLimits = [0.0001 1;...
                             0.0001 1;...
                             -200 200;...
                             0 1];
    fparams.ModelEvidence.integrationMethod = 'quad';
    fparams.ModelEvidence.integrationOptions.dx = 0.1;
    fparams.ModelEvidence.OpenMind = 1;
end