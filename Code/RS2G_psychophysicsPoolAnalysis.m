function [mtp_in, stdtp_in, bias, variance, rmse, WM, WP, B, pval, weber,...
    tsIn, tpIn, trialsIn, ts_in, tp_in, Trials_sorted, estb, ta, G,...
    Llikelihood, BIASs, VARs, lapse, lapseTrials, simbias, simv, RMSEs,...
    LmodelEvidence, notFitLlikelihood, notFitLmodelEvidence,...
    WM_DRIFT, W_INT, ALPHA] = RS2G_psychophysicsPoolAnalysis(d,varargin)
%% RS4GPoolAnalysis
%
%   [mtp_in, stdtp_in, bias, variance, rmse, wm, wp, B, pval, weber, tsIn, tpIn, trialsIn, ts_in, tp_in, Trials_sorted, estb, llikelihood, ta, ta2] = RS4GPoolAnalysis(varargin)   
%
%%

%% Defaults
Fit_default.fittype = 'none';       % No fit by default
Fit_default.trialtypes = 1;
Fit_default.CrossValidation.Type = 'None';
Fit_default.ModelEvidence.method = 'None';
Fit_default.modelUsed = 1;
Fit_default.ObsAct = 0;
bootstrap_default.nbootstraps = NaN;
bootstrap_default.nsamps = NaN;
PlotOpts_default.titles = {'RS1G','RS2G'};
PlotOpts_default.RelativeFigSize = [1/5 1/2 3/5 1/3];
PlotOpts_default.colors = [0 0 1; 1 0 0; 0.6 0.6 0.6; 0 0 0];
TAexpectation_default.method = 'none';

%% Parse input
Parser = inputParser;

addRequired(Parser,'d')     % Data structure
addParameter(Parser,'runs',NaN)     % Runs to analyze
addParameter(Parser,'interval_N',1:2)    % Trial types to analyze
addParameter(Parser,'outlier',Inf)      % Number of standard deviations away from the mean to exclude data as outlier
addParameter(Parser,'trialWin',[1 Inf]) % Trials in each run to analyze ([start# end#])
addParameter(Parser,'Fit',Fit_default)  % Fitting options
addParameter(Parser,'Plot','none')      % Data to plot
addParameter(Parser,'PlotOpts',PlotOpts_default)    % Plotting options
addParameter(Parser,'Save','No')           % Saving options
addParameter(Parser,'WeberFractionCheck',0) % Check fit of weber fraction against data
addParameter(Parser,'Bootstrap',bootstrap_default)  % Bootstrap of bias/variance
addParameter(Parser,'MinMaxTp',[-Inf Inf])      % Minimum and maximum values of t_p to keep
addParameter(Parser,'ConflictType','equal')     % For experiments with cue conflict
addParameter(Parser,'DiffTolerance',2/60)       % Tolerance for difference in sample times before calling it conflict
addParameter(Parser,'TAexpectation',TAexpectation_default)  % For controlling the calculation of the expected value of aim times under a model
addParameter(Parser,'OutlierRejectionRounds',3)

parse(Parser,d,varargin{:})

d = Parser.Results.d;
runs = Parser.Results.runs;
interval_N = Parser.Results.interval_N;
outlier = Parser.Results.outlier;
trialWin = Parser.Results.trialWin;
Fit = Parser.Results.Fit;
Plot = Parser.Results.Plot;
PlotOpts = Parser.Results.PlotOpts;
Save = Parser.Results.Save;
WeberFractionCheck = Parser.Results.WeberFractionCheck;
Bootstrap = Parser.Results.Bootstrap;
MinMaxTp = Parser.Results.MinMaxTp;
ConflictType = Parser.Results.ConflictType;
DiffTolerance = Parser.Results.DiffTolerance;
TAexpectation = Parser.Results.TAexpectation;
OutlierRejectionRounds = Parser.Results.OutlierRejectionRounds;

% Check to see if run information was provided
if isnan(runs)
    runs = 1:d.runs;            % Defaults to pool data from all runs
end

% Set m to interval_Ns
m = interval_N;

% Pull out bootstrap parameters
if any(isnan([Bootstrap.nbootstraps Bootstrap.nsamps]))
    bootflg = 0;
else
    bootflg = 1;
    nbootstraps = Bootstrap.nbootstraps;
    nsamps = Bootstrap.nsamps;
end

% Determine if fitting options have fields for CrossValidation and
% ModelEvdience
if ~isfield(Fit,'CrossValidation')
    Fit.CrossValidation.Type = 'None';
end
if ~isfield(Fit,'ModelEvidence')
    Fit.ModelEvidence.method = 'None';
end
if ~isfield(Fit,'modelUsed')
    Fit.modelUsed = 1;      % Use the first model fit by default
end
if ~isfield(Fit,'ObsAct')
    Fit.ObsAct = 0;
end

%% Analyze data

for i = m
    % Grab the appropriate data
    [ts1{i}, ts2{i}, ts{i}, tp{i}, ~, ~, Trials{i}, correct{i}, ~, ~, ~, N{i}] = RS2G_psychophysics_pooldata(d,'runs',runs,'interval_N',i,'trialWin',trialWin);
    
    % Find data associated with the desired conflict type
    switch ConflictType
        case 'all'
            ts{i} = ts{i};
            Trials{i} = Trials{i};
        case 'equal'
            ts{i} = ts{i}(abs(ts1{i} - ts2{i}) <= DiffTolerance);
            Trials{i} = Trials{i}(abs(ts1{i} - ts2{i}) <= DiffTolerance);
            tp{i} = tp{i}(abs(ts1{i} - ts2{i}) <= DiffTolerance);
        case 'ts1 > ts2'
            ts{i} = ts{i}(ts1{i} > ts2{i} & abs(ts1{i}-ts2{i}) > DiffTolerance);
            Trials{i} = Trials{i}(ts1{i} > ts2{i} & abs(ts1{i}-ts2{i}) > DiffTolerance);
            tp{i} = tp{i}(ts1{i} > ts2{i} & abs(ts1{i}-ts2{i}) > DiffTolerance);
        case 'ts1 < ts2'
            ts{i} = ts{i}(ts1{i} < ts2{i} & abs(ts1{i}-ts2{i}) > DiffTolerance);
            Trials{i} = Trials{i}(ts1{i} < ts2{i} & abs(ts1{i}-ts2{i}) > DiffTolerance);
            tp{i} = tp{i}(ts1{i} < ts2{i} & abs(ts1{i}-ts2{i}) > DiffTolerance);
    end
    tss = unique(ts{i});
    
    % Estimate Weber fraction
    disp('Estimating Weber fraction')
    if isempty(tp{i})
        weber(i) = 0.1;
    elseif sum(tp{i} > 0)
        fun = @(w)(sum(w.^2*tss.^2) - var(tp{i}(tp{i} > 0 & tp{i} >= MinMaxTp(1) & tp{i} <= tss(end)+MinMaxTp(2))));
        weber(i) = lsqnonlin(fun,0.1);
    else
        weber(i) = NaN;
    end
    
    % Sort data by sample time
    errors{i} = [];
    tsIn{i} = [];
    tpIn{i} = [];
    trialsIn{i} = [];
    for ii = 1:length(tss)
        ts_sorted{i}{ii} = ts{i}(ts{i} == tss(ii));
        tp_sorted{i}{ii} = tp{i}(ts{i} == tss(ii));
        Trials_sorted{i}{ii} = Trials{i}(ts{i} == tss(ii));
        
        
        mtp_in(ii,i) = nanmean(tp_sorted{i}{ii}(tp_sorted{i}{ii} >= MinMaxTp(1) & tp_sorted{i}{ii} <= tss(ii)+MinMaxTp(2)));
        
        stdtp_in(ii,i) = nanstd(tp_sorted{i}{ii}(tp_sorted{i}{ii} >= MinMaxTp(1) & tp_sorted{i}{ii} <= tss(ii)+MinMaxTp(2)),1); %weber(i)*tss(ii);
        for jj = 1:OutlierRejectionRounds
            mtp(ii,i) = mtp_in(ii,i);
            
            stdtp(ii,i) = stdtp_in(ii,i);
            
            % Find outliers and recalculate mean and variance
            ins{i}{ii} = find(abs(tp_sorted{i}{ii} - mtp(ii,i)) < outlier*stdtp(ii,i) & tp_sorted{i}{ii} >= MinMaxTp(1) & tp_sorted{i}{ii} <= tss(ii)+MinMaxTp(2));
            
            
            mtp_in(ii,i) = nanmean(tp_sorted{i}{ii}(ins{i}{ii}));
            
            stdtp_in(ii,i) = nanstd(tp_sorted{i}{ii}(ins{i}{ii}),1);
            
        end
        ts_in{i}{ii} = ts_sorted{i}{ii}(ins{i}{ii});
        tp_in{i}{ii} = tp_sorted{i}{ii}(ins{i}{ii});
        trials_in{i}{ii} = Trials_sorted{i}{ii}(ins{i}{ii});
        tsIn{i} = [tsIn{i}; ts_sorted{i}{ii}(ins{i}{ii})];
        tpIn{i} = [tpIn{i}; tp_sorted{i}{ii}(ins{i}{ii})];
        trialsIn{i} = [trialsIn{i}; Trials_sorted{i}{ii}(ins{i}{ii})];
        errors{i} = [errors{i}; tss(ii) - tp_in{i}{ii}];
        
        
        % Weber fraction fit check
        if WeberFractionCheck
            if ~isempty(tp{i})
                figure
                [n bins] = hist(tp_sorted{i}{ii},15);
                bar(bins,n/sum(n))
                hold on
                plot(bins,mvnpdf(bins',mean(tp_sorted{i}{ii}),std(tp_sorted{i}{ii}).^2)*(bins(2)-bins(1)),'r--')
                plot(bins,mvnpdf(bins',mean(tp_sorted{i}{ii}),tss(1).^2*weber(i).^2)*(bins(2)-bins(1)),'r')
                title('Go1')
            end
        end
    end
    
    
    % Test if production times for longest interval are significantly longer than shortest interval
    tpshort = tp_in{i}{1};
    tplong = tp_in{i}{end};
    
    [h, pval(i)] = ttest2(tplong,tpshort,'tail','right');
    
end         % Sets loop

% Fit the data
if ~iscell(Fit.fittype)
    fittype{1} = Fit.fittype;
else
    fittype = Fit.fittype;
end
for fits = 1:length(fittype)
    switch fittype{fits}
        case 'BLSbiased'
            switch Fit.method
                case 'quad'
                    disp('Fitting BLS with bias model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estb'
                                init = [0.1 0.06 NaN];      % Default initial search but with estimate of baseline
                                estb(i) = nanmean(tpIn{i}) - nanmean(tsIn{i});
                                init(3) = estb(i);
                            case 'default'
                                init = [0.1 0.06 0];
                        end
                    end
                    dt = Fit.dx;
                    [WM(:,fits), WP(:,fits,1), B(:,fits), Llikelihood(:,fits)] = BLSbiased_fitter(tsIn(Fit.trialtypes),tpIn(Fit.trialtypes),'InitCond',init,'FitType',Fit.method,[min(tss) max(tss) dt],'N',num2cell(Fit.trialtypes),'ObsAct',Fit.ObsAct);
                    LmodelEvidence(:,fits) = NaN;
                case 'quad_batch'
                    disp('Fitting BLS with bias model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estb'
                                init = [0.1 0.06 NaN];      % Default initial search but with estimate of baseline
                                estb(i) = nanmean(tpIn{i}) - nanmean(tsIn{i});
                                init(3) = estb(i);
                            case 'default'
                                init = [0.1 0.06 0];
                        end
                    end
                    dt = Fit.dx;
                    if isfield(varargin{fitnum+1},'batchsize')
                        batchsize = varargin{fitnum+1}.batchsize;
                    else
                        batchsize = 1000000;
                    end
                    [WM(:,fits), WP(:,fits,1), B(:,fits), Llikelihood(:,fits)] = BLSbiased_fitter(tsIn(Fit.trialtypes),tpIn(Fit.trialtypes),'InitCond',init,'FitType',Fit.method,[min(tss) max(tss) dt batchsize],'N',num2cell(Fit.trialtypes),'ObsAct',Fit.ObsAct);
                    LmodelEvidence(:,fits) = NaN;
                otherwise
                    error('Fitting method not recognized for BLSbiased fitter!')
            end
            G(:,fits) = 1;
            lapse(:,fits) = 0;
            if fits == Fit.modelUsed;
                for i = 1:length(tsIn)
                    lapseTrials{i} = ~ones(size(tsIn{i}));
                end
            end
            WM_DRIFT(:,fits) = WM(:,fits);
            W_INT(:,fits) = WM(:,fits)/sqrt(2);
            ALPHA(:,fits) = ones(size(WM,1),1);
            WP(:,fits,2) = WP(:,fits,1);
            
        case 'gBLS'
            switch Fit.method
                case 'quad'
                    disp('Fitting BLS with gain model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estg'
                                init = [0.1 0.06 NaN];      % Default initial search but with estimate of baseline
                                estg(i) = regress(tpIn{i},tsIn{i});
                                init(3) = estg(i);
                            case 'default'
                                init = [0.1 0.06 1];
                        end
                    end
                    dt = Fit.dx;
                    [WM(:,fits), WP(:,fits,1), G(:,fits), Llikelihood(:,fits)] = gBLS_fitter(tsIn(Fit.trialtypes),tpIn(Fit.trialtypes),'InitCond',init,'FitType',Fit.method,[min(tss) max(tss) dt],'N',num2cell(Fit.trialtypes),'ObsAct',Fit.ObsAct);
                    LmodelEvidence(:,fits) = NaN;
                case 'quad_batch'
                    disp('Fitting BLS with gain model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estg'
                                init = [0.1 0.06 NaN];      % Default initial search but with estimate of baseline
                                estg(i) = regress(tpIn{i},tsIn{i});
                                init(3) = estg(i);
                            case 'default'
                                init = [0.1 0.06 1];
                        end
                    end
                    dt = Fit.dx;
                    if isfield(varargin{fitnum+1},'batchsize')
                        batchsize = varargin{fitnum+1}.batchsize;
                    else
                        batchsize = 1000000;
                    end
                    [WM(:,fits), WP(:,fits,1), G(:,fits), Llikelihood(:,fits)] = gBLS_fitter(tsIn(Fit.trialtypes),tpIn(Fit.trialtypes),'InitCond',init,'FitType',Fit.method,[min(tss) max(tss) dt batchsize],'N',num2cell(Fit.trialtypes),'ObsAct',Fit.ObsAct);
                    LmodelEvidence = NaN;
                otherwise
                    error('Fitting method not recognized for gBLS fitter!')
            end
            B(:,fits) = 0;
            lapse(:,fits) = 0;
            if fits == Fit.modelUsed;
                for i = 1:length(tsIn)
                    lapseTrials{i} = ~ones(size(tsIn{i}));
                end
            end
            WM_DRIFT(:,fits) = WM(:,fits);
            W_INT(:,fits) = WM(:,fits)/sqrt(2);
            ALPHA(:,fits) = ones(size(WM,1),1);
            WP(:,fits,2) = WP(:,fits,1);
            
        case 'BLSbiasedLapse'
            switch Fit.method
                case 'quad'
                    % Fit the BLS model with bias and lapses to the data
                    disp('Fitting BLS with bias and lapses model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estb'
                                init = [0.1 0.06 NaN 0.05];      % Default initial search but with estimate of baseline
                                estb(i) = nanmean(tpIn{i}) - nanmean(tsIn{i});
                                init(3) = estb(i);
                            case 'default'
                                init = [0.1 0.06 0 0.05];
                        end
                    end
                    dt = Fit.dx;
                    LapseSupport = [MinMaxTp(1) max(tss)+MinMaxTp(2)];
                    [WM(:,fits), WP(:,fits,1), B(:,fits), lapse(:,fits), Llikelihood(:,fits), LmodelEvidence(:,fits)] = BLSbiasedLapse_fitter(tsIn(Fit.trialtypes),tpIn(Fit.trialtypes),'InitCond',init,'FitType',Fit.method,[min(tss) max(tss) dt],'N',num2cell(Fit.trialtypes),'LapseSupport',LapseSupport,...
                        'CrossValidation',Fit.CrossValidation,'ModelEvidence',Fit.ModelEvidence,'ObsAct',Fit.ObsAct);
                    
                    % Find likelihood of model on left out condition, if it exists
                    for i = m
                        fitted(i) = any(i == Fit.trialtypes);
                    end
                    if any(~fitted)
                        [notFitLlikelihood(:,fits), notFitLmodelEvidence(:,fits)] = BLSbiasedLapse_Validator(tsIn(~fitted),tpIn(~fitted),mean(WM(:,fits)),mean(WP(:,fits,1)),mean(B(:,fits)),mean(lapse(:,fits)),...
                            'N',num2cell(m(~fitted)),'LapseSupport',LapseSupport,'ModelEvidence',Fit.ModelEvidence,'FitType',Fit.method,[min(tss) max(tss) dt],'ObsAct',Fit.ObsAct); 
                    else
                        notFitLlikelihood(:,fits) = NaN;
                        notFitLmodelEvidence(:,fits) = NaN;
                    end
                    
                    % Identify lapse trials
                    if fits == Fit.modelUsed
                        estimator.type = 'BLS';
                        estimator.ObsAct = Fit.ObsAct;
                        estimator.wy = mean(WM(:,fits));
                        for i = m
                            [~, ~, loglike, ~, like] = prob_tp_take_ts_wm_wp(tpIn{i}-mean(B(:,fits)),tsIn{i},mean(WM(:,fits)),mean(WP(:,fits,1)),i,'estimator',estimator);
                            loglikeLapse = log(mean(lapse(:,fits))/(LapseSupport(2)-LapseSupport(1)));
                            lapseTrials{i} = log((1-mean(lapse(:,fits)))*like) < loglikeLapse;
                        end
                    end
                    
                    % Recalculate mtp and stdtp
                    for i = m
                        for ii = 1:length(tss)
                            mtp_in(ii,i) = mean(tpIn{i}(tsIn{i} == tss(ii) & ~lapseTrials{i}));
                            stdtp_in(ii,i) = std(tpIn{i}(tsIn{i} == tss(ii) & ~lapseTrials{i}));
                        end
                        errors{i} = tsIn{i}(~lapseTrials{i}) - tpIn{i}(~lapseTrials{i});
                    end
                    
                case 'quad_batch'
                    disp('Fitting BLS with bias model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estb'
                                init = [0.1 0.06 NaN];      % Default initial search but with estimate of baseline
                                estb(i) = nanmean(tpIn{i}) - nanmean(tsIn{i});
                                init(3) = estb(i);
                            case 'default'
                                init = [0.1 0.06 0];
                        end
                    end
                    dt = Fit.dx;
                    if isfield(varargin{fitnum+1},'batchsize')
                        batchsize = varargin{fitnum+1}.batchsize;
                    else
                        batchsize = 1000000;
                    end
                    LapseSupport = [MinMaxTp(1) max(tss)+MinMaxTp(2)];
                    [WM(:,fits), WP(:,fits,1), B(:,fits), lapse(:,fits), Llikelihood(:,fits), LmodelEvidence(:,fits)] = BLSbiasedLapse_fitter(tsIn(Fit.trialtypes),tpIn(Fit.trialtypes),'InitCond',init,'FitType',Fit.method,[min(tss) max(tss) dt batchsize],'N',num2cell(Fit.trialtypes),'LapseSupport',LapseSupport,...
                        'CrossValidation',Fit.CrossValidation,'ModelEvidence',Fit.ModelEvidence,'ObsAct',Fit.ObsAct);
                    
                    % Find likelihood of model on left out condition, if it exists
                    for i = m
                        fitted(i) = any(i == Fit.trialtypes);
                    end
                    if any(~fitted)
                        [notFitLlikelihood(:,fits), notFitLmodelEvidence(:,fits)] = BLSbiasedLapse_Validator(tsIn(~fitted),tpIn(~fitted),mean(WM(:,fits)),mean(WP(:,fits,1)),mean(B(:,fits)),mean(lapse(:,fits)),...
                            'N',num2cell(m(~fitted)),'LapseSupport',LapseSupport,'ModelEvidence',Fit.ModelEvidence,'FitType',Fit.method,[min(tss) max(tss) dt],'ObsAct',Fit.ObsAct); 
                    else
                        notFitLlikelihood(:,fits) = NaN;
                        notFitLmodelEvidence(:,fits) = NaN;
                    end
                    
                    % Identify lapse trials
                    if fits == Fit.modelUsed
                        estimator.type = 'BLS';
                        estimator.ObsAct = Fit.ObsAct;
                        estimator.wy = mean(WM(:,fits));
                        for i = m
                            [~, ~, loglike, ~, like] = prob_tp_take_ts_wm_wp(tpIn{i}-mean(B(:,fits)),tsIn{i},mean(WM(:,fits)),mean(WP(:,fits,1)),i,'estimator',estimator);
                            loglikeLapse = log(mean(lapse(:,fits))/(LapseSupport(2)-LapseSupport(1)));
                            lapseTrials{i} = log((1-mean(lapse(:,fits)))*like) < loglikeLapse;
                        end
                    end
                    
                    % Recalculate mtp and stdtp
                    for i = m
                        for ii = 1:length(tss)
                            mtp_in(ii,i) = mean(tpIn{i}(tsIn{i} == tss(ii) & ~lapseTrials{i}));
                            stdtp_in(ii,i) = std(tpIn{i}(tsIn{i} == tss(ii) & ~lapseTrials{i}));
                        end
                        errors{i} = tsIn{i}(~lapseTrials{i}) - tpIn{i}(~lapseTrials{i});
                    end
                    
                otherwise
                    error('Fitting method not recognized for BLSbiased fitter!')
            end
            G(:,fits) = 1;
            WM_DRIFT(:,fits) = WM(:,fits);
            W_INT(:,fits) = WM(:,fits)/sqrt(2);
            ALPHA(:,fits) = ones(size(WM,1),1);
            WP(:,fits,2) = WP(:,fits,1);
        
            
        case 'MAPbiasedLapse'
            switch Fit.method
                case 'quad'
                    % Fit the BLS model with bias and lapses to the data
                    disp('Fitting MAP with bias and lapses model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estb'
                                init = [0.1 0.06 NaN 0.05];      % Default initial search but with estimate of baseline
                                estb(i) = nanmean(tpIn{i}) - nanmean(tsIn{i});
                                init(3) = estb(i);
                            case 'default'
                                init = [0.1 0.06 0 0.05];
                        end
                    end
                    dt = Fit.dx;
                    LapseSupport = [MinMaxTp(1) max(tss)+MinMaxTp(2)];
                    [WM(:,fits), WP(:,fits,1), B(:,fits), lapse(:,fits), Llikelihood(:,fits), LmodelEvidence(:,fits)] = MAPbiasedLapse_fitter(tsIn(Fit.trialtypes),tpIn(Fit.trialtypes),'InitCond',init,'FitType',Fit.method,[min(tss) max(tss) dt],'N',num2cell(Fit.trialtypes),'LapseSupport',LapseSupport,...
                        'CrossValidation',Fit.CrossValidation,'ModelEvidence',Fit.ModelEvidence);
                    
                    % Find likelihood of model on left out condition, if it exists
                    for i = m
                        fitted(i) = any(i == Fit.trialtypes);
                    end
                    if any(~fitted)
                        [notFitLlikelihood(:,fits), notFitLmodelEvidence(:,fits)] = MAPbiasedLapse_Validator(tsIn(~fitted),tpIn(~fitted),mean(WM(:,fits)),mean(WP(:,fits,1)),mean(B(:,fits)),mean(lapse(:,fits)),...
                            'N',num2cell(m(~fitted)),'LapseSupport',LapseSupport,'ModelEvidence',Fit.ModelEvidence,'FitType',Fit.method,[min(tss) max(tss) dt]); 
                    else
                        notFitLlikelihood(:,fits) = NaN;
                        notFitLmodelEvidence(:,fits) = NaN;
                    end
                    
                    % Identify lapse trials
                    if fits == Fit.modelUsed
                        for i = m
                            [~, ~, loglike, ~, like] = prob_tp_take_ts_wm_wp(tpIn{i}-mean(B(:,fits)),tsIn{i},mean(WM(:,fits)),mean(WP(:,fits,1)),i);
                            loglikeLapse = log(mean(lapse(:,fits))/(LapseSupport(2)-LapseSupport(1)));
                            lapseTrials{i} = log((1-mean(lapse(:,fits)))*like) < loglikeLapse;
                        end
                    end
                    
                    % Recalculate mtp and stdtp
                    for i = m
                        for ii = 1:length(tss)
                            mtp_in(ii,i) = mean(tpIn{i}(tsIn{i} == tss(ii) & ~lapseTrials{i}));
                            stdtp_in(ii,i) = std(tpIn{i}(tsIn{i} == tss(ii) & ~lapseTrials{i}));
                        end
                        errors{i} = tsIn{i}(~lapseTrials{i}) - tpIn{i}(~lapseTrials{i});
                    end
                    
                case 'quad_batch'
                    disp('Fitting MAP with bias model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estb'
                                init = [0.1 0.06 NaN];      % Default initial search but with estimate of baseline
                                estb(i) = nanmean(tpIn{i}) - nanmean(tsIn{i});
                                init(3) = estb(i);
                            case 'default'
                                init = [0.1 0.06 0];
                        end
                    end
                    dt = Fit.dx;
                    if isfield(varargin{fitnum+1},'batchsize')
                        batchsize = varargin{fitnum+1}.batchsize;
                    else
                        batchsize = 1000000;
                    end
                    LapseSupport = [MinMaxTp(1) max(tss)+MinMaxTp(2)];
                    [WM(:,fits), WP(:,fits,1), B(:,fits), lapse(:,fits), Llikelihood(:,fits), LmodelEvidence(:,fits)] = MAPbiasedLapse_fitter(tsIn(Fit.trialtypes),tpIn(Fit.trialtypes),'InitCond',init,'FitType',Fit.method,[min(tss) max(tss) dt batchsize],'N',num2cell(Fit.trialtypes),'LapseSupport',LapseSupport,...
                        'CrossValidation',Fit.CrossValidation,'ModelEvidence',Fit.ModelEvidence);
                    
                    % Find likelihood of model on left out condition, if it exists
                    for i = m
                        fitted(i) = any(i == Fit.trialtypes);
                    end
                    if any(~fitted)
                        [notFitLlikelihood(:,fits), notFitLmodelEvidence(:,fits)] = BLSbiasedLapse_Validator(tsIn(~fitted),tpIn(~fitted),mean(WM(:,fits)),mean(WP(:,fits,1)),mean(B(:,fits)),mean(lapse(:,fits)),...
                            'N',num2cell(m(~fitted)),'LapseSupport',LapseSupport,'ModelEvidence',Fit.ModelEvidence,'FitType',Fit.method,[min(tss) max(tss) dt]); 
                    else
                        notFitLlikelihood(:,fits) = NaN;
                        notFitLmodelEvidence(:,fits) = NaN;
                    end
                    
                    % Identify lapse trials
                    if fits == Fit.modelUsed;
                        for i = m
                            [~, ~, loglike, ~, like] = prob_tp_take_ts_wm_wp(tpIn{i}-mean(B(:,fits)),tsIn{i},mean(WM(:,fits)),mean(WP(:,fits,1)),i);
                            loglikeLapse = log(mean(lapse(:,fits))/(LapseSupport(2)-LapseSupport(1)));
                            lapseTrials{i} = log((1-mean(lapse(:,fits)))*like) < loglikeLapse;
                        end
                    end
                    
                    % Recalculate mtp and stdtp
                    for i = m
                        for ii = 1:length(tss)
                            mtp_in(ii,i) = mean(tpIn{i}(tsIn{i} == tss(ii) & ~lapseTrials{i}));
                            stdtp_in(ii,i) = std(tpIn{i}(tsIn{i} == tss(ii) & ~lapseTrials{i}));
                        end
                        errors{i} = tsIn{i}(~lapseTrials{i}) - tpIn{i}(~lapseTrials{i});
                    end
                    
                otherwise
                    error('Fitting method not recognized for BLSbiased fitter!')
            end
            G(:,fits) = 1;
            WM_DRIFT(:,fits) = WM(:,fits);
            W_INT(:,fits) = WM(:,fits)/sqrt(2);
            ALPHA(:,fits) = ones(size(WM,1),1);
            WP(:,fits,2) = WP(:,fits,1);
            
        case 'ObsActBiasedLapse'
            switch Fit.method
                case 'quad'
                    % Fit the BLS model with bias and lapses to the data
                    disp('Fitting BLS with bias and lapses model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estb'
                                init = [0.1 0.06 NaN 0.05];      % Default initial search but with estimate of baseline
                                estb(i) = nanmean(tpIn{i}) - nanmean(tsIn{i});
                                init(3) = estb(i);
                            case 'default'
                                init = [0.1 0.06 0 0.05];
                        end
                    end
                    dt = Fit.dx;
                    LapseSupport = [MinMaxTp(1) max(tss)+MinMaxTp(2)];
                    [WM(:,fits), WP(:,fits,1), B(:,fits), lapse(:,fits), Llikelihood(:,fits), LmodelEvidence(:,fits)] = ObsActBiasedLapse_fitter(tsIn(Fit.trialtypes),tpIn(Fit.trialtypes),'InitCond',init,'FitType',Fit.method,[min(tss) max(tss) dt],'N',num2cell(Fit.trialtypes),'LapseSupport',LapseSupport,...
                        'CrossValidation',Fit.CrossValidation,'ModelEvidence',Fit.ModelEvidence);
                    
                    % Find likelihood of model on left out condition, if it exists
                    for i = m
                        fitted(i) = any(i == Fit.trialtypes);
                    end
                    if any(~fitted)
                        [notFitLlikelihood(:,fits), notFitLmodelEvidence(:,fits)] = ObsActBiasedLapse_Validator(tsIn(~fitted),tpIn(~fitted),mean(WM(:,fits)),mean(WP(:,fits,1)),mean(B(:,fits)),mean(lapse(:,fits)),...
                            'N',num2cell(m(~fitted)),'LapseSupport',LapseSupport,'ModelEvidence',Fit.ModelEvidence,'FitType',Fit.method,[min(tss) max(tss) dt]); 
                    else
                        notFitLlikelihood(:,fits) = NaN;
                        notFitLmodelEvidence(:,fits) = NaN;
                    end
                    
                    % Identify lapse trials
                    if fits == Fit.modelUsed
                        for i = m
                            [~, ~, loglike, ~, like] = prob_tp_take_ts_wm_wp(tpIn{i}-mean(B(:,fits)),tsIn{i},mean(WM(:,fits)),mean(WP(:,fits,1)),i);
                            loglikeLapse = log(mean(lapse(:,fits))/(LapseSupport(2)-LapseSupport(1)));
                            lapseTrials{i} = log((1-mean(lapse(:,fits)))*like) < loglikeLapse;
                        end
                    end
                    
                    % Recalculate mtp and stdtp
                    for i = m
                        for ii = 1:length(tss)
                            mtp_in(ii,i) = mean(tpIn{i}(tsIn{i} == tss(ii) & ~lapseTrials{i}));
                            stdtp_in(ii,i) = std(tpIn{i}(tsIn{i} == tss(ii) & ~lapseTrials{i}));
                        end
                        errors{i} = tsIn{i}(~lapseTrials{i}) - tpIn{i}(~lapseTrials{i});
                    end
                    
                case 'quad_batch'
                    disp('Fitting BLS with bias model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estb'
                                init = [0.1 0.06 NaN];      % Default initial search but with estimate of baseline
                                estb(i) = nanmean(tpIn{i}) - nanmean(tsIn{i});
                                init(3) = estb(i);
                            case 'default'
                                init = [0.1 0.06 0];
                        end
                    end
                    dt = Fit.dx;
                    if isfield(varargin{fitnum+1},'batchsize')
                        batchsize = varargin{fitnum+1}.batchsize;
                    else
                        batchsize = 1000000;
                    end
                    LapseSupport = [MinMaxTp(1) max(tss)+MinMaxTp(2)];
                    [WM(:,fits), WP(:,fits,1), B(:,fits), lapse(:,fits), Llikelihood(:,fits), LmodelEvidence(:,fits)] = ObsActBiasedLapse_fitter(tsIn(Fit.trialtypes),tpIn(Fit.trialtypes),'InitCond',init,'FitType',Fit.method,[min(tss) max(tss) dt batchsize],'N',num2cell(Fit.trialtypes),'LapseSupport',LapseSupport,...
                        'CrossValidation',Fit.CrossValidation,'ModelEvidence',Fit.ModelEvidence);
                    
                    % Find likelihood of model on left out condition, if it exists
                    for i = m
                        fitted(i) = any(i == Fit.trialtypes);
                    end
                    if any(~fitted)
                        [notFitLlikelihood(:,fits), notFitLmodelEvidence(:,fits)] = ObsActBiasedLapse_Validator(tsIn(~fitted),tpIn(~fitted),mean(WM(:,fits)),mean(WP(:,fits,1)),mean(B(:,fits)),mean(lapse(:,fits)),...
                            'N',num2cell(m(~fitted)),'LapseSupport',LapseSupport,'ModelEvidence',Fit.ModelEvidence,'FitType',Fit.method,[min(tss) max(tss) dt]); 
                    else
                        notFitLlikelihood(:,fits) = NaN;
                        notFitLmodelEvidence(:,fits) = NaN;
                    end
                    
                    % Identify lapse trials
                    if fits == Fit.modelUsed
                        for i = m
                            [~, ~, loglike, ~, like] = prob_tp_take_ts_wm_wp(tpIn{i}-mean(B(:,fits)),tsIn{i},mean(WM(:,fits)),mean(WP(:,fits,1)),i);
                            loglikeLapse = log(mean(lapse(:,fits))/(LapseSupport(2)-LapseSupport(1)));
                            lapseTrials{i} = log((1-mean(lapse(:,fits)))*like) < loglikeLapse;
                        end
                    end
                    
                    % Recalculate mtp and stdtp
                    for i = m
                        for ii = 1:length(tss)
                            mtp_in(ii,i) = mean(tpIn{i}(tsIn{i} == tss(ii) & ~lapseTrials{i}));
                            stdtp_in(ii,i) = std(tpIn{i}(tsIn{i} == tss(ii) & ~lapseTrials{i}));
                        end
                        errors{i} = tsIn{i}(~lapseTrials{i}) - tpIn{i}(~lapseTrials{i});
                    end
                    
                otherwise
                    error('Fitting method not recognized for BLSbiased fitter!')
            end
            G(:,fits) = 1;
            WM_DRIFT(:,fits) = WM(:,fits);
            W_INT(:,fits) = WM(:,fits)/sqrt(2);
            ALPHA(:,fits) = ones(size(WM,1),1);
            WP(:,fits,2) = WP(:,fits,1);
            
        case 'aveMeasurements'
            switch Fit.method
                case 'quad'
                    % Fit the aveMeas model with bias and lapses to the data
                    disp('Fitting measurement averaging model with bias and lapses model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estb'
                                init = [0.1 0.06 NaN 0.05];      % Default initial search but with estimate of baseline
                                estb(i) = nanmean(tpIn{i}) - nanmean(tsIn{i});
                                init(3) = estb(i);
                            case 'default'
                                init = [0.1 0.06 0 0.05];
                        end
                    end
                    dt = Fit.dx;
                    LapseSupport = [MinMaxTp(1) max(tss)+MinMaxTp(2)];
                    [WM(:,fits), WP(:,fits,1), B(:,fits), lapse(:,fits), Llikelihood(:,fits), LmodelEvidence(:,fits)] = aveMeasbiasedLapse_fitter(tsIn(Fit.trialtypes),tpIn(Fit.trialtypes),'InitCond',init,'FitType',Fit.method,[min(tss) max(tss) dt],'N',num2cell(Fit.trialtypes),'LapseSupport',LapseSupport,...
                        'CrossValidation',Fit.CrossValidation,'ModelEvidence',Fit.ModelEvidence,'ObsAct',Fit.ObsAct);
                    
                    % Find likelihood of model on left out condition, if it exists
                    for i = m
                        fitted(i) = any(i == Fit.trialtypes);
                    end
                    if any(~fitted)
                        [notFitLlikelihood(:,fits), notFitLmodelEvidence(:,fits)] = aveMeasbiasedLapse_Validator(tsIn(~fitted),tpIn(~fitted),mean(WM(:,fits)),mean(WP(:,fits,1)),mean(B(:,fits)),mean(lapse(:,fits)),...
                            'N',num2cell(m(~fitted)),'LapseSupport',LapseSupport,'ModelEvidence',Fit.ModelEvidence,'FitType',Fit.method,[min(tss) max(tss) dt],'ObsAct',Fit.ObsAct); 
                    else
                        notFitLlikelihood(:,fits) = NaN;
                        notFitLmodelEvidence(:,fits) = NaN;
                    end
                    
                    % Identify lapse trials
                    if fits == Fit.modelUsed
                        for i = m
                            estimator.type = 'weightedMean';
                            estimator.weights = ones(1,Fit.trialtypes(i))/Fit.trialtypes(i);
                            estimator.ObsAct = Fit.ObsAct;
                            estimator.wy = mean(WM(:,fits));
                            [~, ~, loglike, ~, like] = prob_tp_take_ts_wm_wp(tpIn{i}-mean(B(:,fits)),tsIn{i},mean(WM(:,fits)),mean(WP(:,fits,1)),i,'estimator',estimator);
                            loglikeLapse = log(mean(lapse(:,fits))/(LapseSupport(2)-LapseSupport(1)));
                            lapseTrials{i} = log((1-mean(lapse(:,fits)))*like) < loglikeLapse;
                        end
                    end
                    
                    % Recalculate mtp and stdtp
                    for i = m
                        for ii = 1:length(tss)
                            mtp_in(ii,i) = mean(tpIn{i}(tsIn{i} == tss(ii) & ~lapseTrials{i}));
                            stdtp_in(ii,i) = std(tpIn{i}(tsIn{i} == tss(ii) & ~lapseTrials{i}));
                        end
                        errors{i} = tsIn{i}(~lapseTrials{i}) - tpIn{i}(~lapseTrials{i});
                    end
                    
                case 'quad_batch'
                    disp('Fitting measurement averaging model with bias model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estb'
                                init = [0.1 0.06 NaN];      % Default initial search but with estimate of baseline
                                estb(i) = nanmean(tpIn{i}) - nanmean(tsIn{i});
                                init(3) = estb(i);
                            case 'default'
                                init = [0.1 0.06 0];
                        end
                    end
                    dt = Fit.dx;
                    if isfield(varargin{fitnum+1},'batchsize')
                        batchsize = varargin{fitnum+1}.batchsize;
                    else
                        batchsize = 1000000;
                    end
                    LapseSupport = [MinMaxTp(1) max(tss)+MinMaxTp(2)];
                    [WM(:,fits), WP(:,fits,1), B(:,fits), lapse(:,fits), Llikelihood(:,fits), LmodelEvidence(:,fits)] = aveMeasbiasedLapse_fitter(tsIn(Fit.trialtypes),tpIn(Fit.trialtypes),'InitCond',init,'FitType',Fit.method,[min(tss) max(tss) dt batchsize],'N',num2cell(Fit.trialtypes),'LapseSupport',LapseSupport,...
                        'CrossValidation',Fit.CrossValidation,'ModelEvidence',Fit.ModelEvidence,'ObsAct',Fit.ObsAct);
                    
                    % Find likelihood of model on left out condition, if it exists
                    for i = m
                        fitted(i) = any(i == Fit.trialtypes);
                    end
                    if any(~fitted)
                        [notFitLlikelihood(:,fits), notFitLmodelEvidence(:,fits)] = aveMeasbiasedLapse_Validator(tsIn(~fitted),tpIn(~fitted),mean(WM(:,fits)),mean(WP(:,fits,1)),mean(B(:,fits)),mean(lapse(:,fits)),...
                            'N',num2cell(m(~fitted)),'LapseSupport',LapseSupport,'ModelEvidence',Fit.ModelEvidence,'FitType',Fit.method,[min(tss) max(tss) dt],'ObsAct',Fit.ObsAct); 
                    else
                        notFitLlikelihood(:,fits) = NaN;
                        notFitLmodelEvidence(:,fits) = NaN;
                    end
                    
                    % Identify lapse trials
                    if fits == Fit.modelUsed
                        for i = m
                            estimator.type = 'weightedMean';
                            estimator.weights = ones(1,Fit.trialtypes(i))/Fit.trialtypes(i);
                            estimator.wy = mean(WM(:,fits));
                            [~, ~, loglike, ~, like] = prob_tp_take_ts_wm_wp(tpIn{i}-mean(B(:,fits)),tsIn{i},mean(WM(:,fits)),mean(WP(:,fits,1)),i,'estimator',estimator);
                            loglikeLapse = log(mean(lapse(:,fits))/(LapseSupport(2)-LapseSupport(1)));
                            lapseTrials{i} = log((1-mean(lapse(:,fits)))*like) < loglikeLapse;
                        end
                    end
                    
                    % Recalculate mtp and stdtp
                    for i = m
                        for ii = 1:length(tss)
                            mtp_in(ii,i) = mean(tpIn{i}(tsIn{i} == tss(ii) & ~lapseTrials{i}));
                            stdtp_in(ii,i) = std(tpIn{i}(tsIn{i} == tss(ii) & ~lapseTrials{i}));
                        end
                        errors{i} = tsIn{i}(~lapseTrials{i}) - tpIn{i}(~lapseTrials{i});
                    end
                    
                otherwise
                    error('Fitting method not recognized for BLSbiased fitter!')
            end
            G(:,fits) = 1;
            WM_DRIFT(:,fits) = WM(:,fits);
            W_INT(:,fits) = WM(:,fits)/sqrt(2);
            ALPHA(:,fits) = ones(size(WM,1),1);
            WP(:,fits,2) = WP(:,fits,1);
        
        case 'SubOptMemBias'
            switch Fit.method
                case 'quad'
                    % Fit the aveMeas model with bias and lapses to the data
                    disp('Fitting measurement SubOptMemBias model with bias and lapses model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estb'
                                init = [0.1 0.06 0.1 0.1/sqrt(2) NaN 0.05];      % Default initial search but with estimate of baseline
                                estb(i) = nanmean(tpIn{i}) - nanmean(tsIn{i});
                                init(5) = estb(i);
                            case 'default'
                                init = [0.1 0.06 0.1 0.1/sqrt(2) 0 0.05];
                        end
                    end
                    dt = Fit.dx;
                    LapseSupport = [MinMaxTp(1) max(tss)+MinMaxTp(2)];
                    [WM(:,fits), WP(:,fits,1), WM_DRIFT(:,fits), W_INT(:,fits),...
                     B(:,fits), lapse(:,fits),Llikelihood(:,fits),...
                     LmodelEvidence(:,fits)] = SubOptMemBiasbiasedLapse_fitter(...
                        tsIn(Fit.trialtypes),tpIn(Fit.trialtypes),...
                        'InitCond',init,'FitType',Fit.method,[min(tss) max(tss) dt],...
                        'N',num2cell(Fit.trialtypes),'LapseSupport',LapseSupport,...
                        'CrossValidation',Fit.CrossValidation,...
                        'ModelEvidence',Fit.ModelEvidence,'ObsAct',Fit.ObsAct,...
                        'InequalityBound',true);
                    
                    % Find likelihood of model on left out condition, if it exists
                    for i = m
                        fitted(i) = any(i == Fit.trialtypes);
                    end
                    if any(~fitted)
                        [notFitLlikelihood(:,fits), notFitLmodelEvidence(:,fits)] = ...
                         SubOptMemBiasbiasedLapse_Validator(...
                            tsIn(~fitted),tpIn(~fitted),...
                            mean(WM(:,fits)),mean(WP(:,fits,1)),...
                            mean(WM_DRIFT(:,fits)),mean(W_INT(:,fits)),...
                            mean(B(:,fits)),mean(lapse(:,fits)),...
                            'N',num2cell(m(~fitted)),...
                            'LapseSupport',LapseSupport,...
                            'ModelEvidence',Fit.ModelEvidence,...
                            'FitType',Fit.method,[min(tss) max(tss) dt],'ObsAct',Fit.ObsAct);
                    else
                        notFitLlikelihood(:,fits) = NaN;
                        notFitLmodelEvidence(:,fits) = NaN;
                    end
                    
                    % Identify lapse trials
                    if fits == Fit.modelUsed
                        for i = m
                            estimator.type = 'SubOptMemBias';
                            estimator.wm_drift = mean(WM_DRIFT(:,fits));
                            estimator.w_int = mean(W_INT(:,fits));
                            estimator.ObsAct = Fit.ObsAct;
                            estimator.wy = mean(WM(:,fits));
                            [~, ~, loglike, ~, like] = prob_tp_take_ts_wm_wp(...
                                tpIn{i}-mean(B(:,fits)),tsIn{i},...
                                mean(WM(:,fits)),mean(WP(:,fits,1)),i,...
                                'estimator',estimator);
                            loglikeLapse = log(mean(lapse(:,fits))/(LapseSupport(2)-LapseSupport(1)));
                            lapseTrials{i} = log((1-mean(lapse(:,fits)))*like) < loglikeLapse;
                        end
                    end
                    
                    % Recalculate mtp and stdtp
                    for i = m
                        for ii = 1:length(tss)
                            mtp_in(ii,i) = mean(tpIn{i}(tsIn{i} == tss(ii) & ~lapseTrials{i}));
                            stdtp_in(ii,i) = std(tpIn{i}(tsIn{i} == tss(ii) & ~lapseTrials{i}));
                        end
                        errors{i} = tsIn{i}(~lapseTrials{i}) - tpIn{i}(~lapseTrials{i});
                    end
                    
                case 'quad_batch'
                    error('Quad_batch not yet supported!')
                    
                otherwise
                    error('Fitting method not recognized for BLSbiased fitter!')
            end
            G(:,fits) = 1;
            ALPHA(:,fits) = ones(size(WM,1),1);
            WP(:,fits,2) = WP(:,fits,1);
            
        case 'NestedModel'
            switch Fit.method
                case 'quad'
                    % Fit the aveMeas model with bias and lapses to the data
                    disp('Fitting NestedModel with bias and lapses model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estb'
                                init = [0.1 0.06 0.1 0.1/sqrt(2) NaN 0.05];      % Default initial search but with estimate of baseline
                                estb(i) = nanmean(tpIn{i}) - nanmean(tsIn{i});
                                init(5) = estb(i);
                            case 'default'
                                init = [0.1 0.06 0.1 0.1/sqrt(2) 0 0.05];
                        end
                    end
                    dt = Fit.dx;
                    LapseSupport = [MinMaxTp(1) max(tss)+MinMaxTp(2)];
                    [WM(:,fits), WP(:,fits,1), WM_DRIFT(:,fits), W_INT(:,fits),...
                        B(:,fits), lapse(:,fits), ALPHA(:,fits), Llikelihood(:,fits),...
                        LmodelEvidence(:,fits)] = NestedModelbiasedLapse_fitter(...
                        tsIn(Fit.trialtypes),tpIn(Fit.trialtypes),...
                        'InitCond',init,'FitType',Fit.method,[min(tss) max(tss) dt],...
                        'N',num2cell(Fit.trialtypes),'LapseSupport',LapseSupport,...
                        'CrossValidation',Fit.CrossValidation,...
                        'ModelEvidence',Fit.ModelEvidence,'ObsAct',Fit.ObsAct);
                    
                    % Find likelihood of model on left out condition, if it exists
                    for i = m
                        fitted(i) = any(i == Fit.trialtypes);
                    end
                    if any(~fitted)
                        [notFitLlikelihood(:,fits), notFitLmodelEvidence(:,fits)] = ...
                            SubOptMemBiasbiasedLapse_Validator(...
                            tsIn(~fitted),tpIn(~fitted),...
                            mean(WM(:,fits)),mean(WP(:,fits,1)),...
                            mean(WM_DRIFT(:,fits)),mean(W_INT(:,fits)),...
                            mean(B(:,fits)),mean(lapse(:,fits)),...
                            'N',num2cell(m(~fitted)),...
                            'LapseSupport',LapseSupport,...
                            'ModelEvidence',Fit.ModelEvidence,...
                            'FitType',Fit.method,[min(tss) max(tss) dt],'ObsAct',Fit.ObsAct);
                    else
                        notFitLlikelihood(:,fits) = NaN;
                        notFitLmodelEvidence(:,fits) = NaN;
                    end
                    
                    % Identify lapse trials
                    if fits == Fit.modelUsed
                        for i = m
                            estimator.type = 'SubOptMemBias';
                            estimator.wm_drift = mean(WM_DRIFT(:,fits));
                            estimator.w_int = mean(W_INT(:,fits));
                            estimator.ObsAct = Fit.ObsAct;
                            estimator.wy = mean(WM(:,fits));
                            [~, ~, loglike, ~, like] = prob_tp_take_ts_wm_wp(...
                                tpIn{i}-mean(B(:,fits)),tsIn{i},...
                                mean(WM(:,fits)),mean(WP(:,fits,1)),i,...
                                'estimator',estimator);
                            loglikeLapse = log(mean(lapse(:,fits))/(LapseSupport(2)-LapseSupport(1)));
                            lapseTrials{i} = log((1-mean(lapse(:,fits)))*like) < loglikeLapse;
                        end
                    end
                    
                    % Recalculate mtp and stdtp
                    for i = m
                        for ii = 1:length(tss)
                            mtp_in(ii,i) = mean(tpIn{i}(tsIn{i} == tss(ii) & ~lapseTrials{i}));
                            stdtp_in(ii,i) = std(tpIn{i}(tsIn{i} == tss(ii) & ~lapseTrials{i}));
                        end
                        errors{i} = tsIn{i}(~lapseTrials{i}) - tpIn{i}(~lapseTrials{i});
                    end
                    
                case 'quad_batch'
                    error('Quad_batch not yet supported!')
                    
                otherwise
                    error('Fitting method not recognized for BLSbiased fitter!')
            end
            G(:,fits) = 1;
            WP(:,fits,2) = WP(:,fits,1);
            
            
        case 'WeightedLinear'
            switch Fit.method
                case 'quad'
                    % Fit the aveMeas model with bias and lapses to the data
                    disp('Fitting measurement averaging model with bias and lapses model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estb'
                                init = [0.1 0.06 NaN 0.05];      % Default initial search but with estimate of baseline
                                estb(i) = nanmean(tpIn{i}) - nanmean(tsIn{i});
                                init(3) = estb(i);
                            case 'default'
                                init = [0.1 0.06 0 0.05];
                        end
                    end
                    dt = Fit.dx;
                    LapseSupport = [MinMaxTp(1) max(tss)+MinMaxTp(2)];
                    [WM(:,fits), WP(:,fits,1), B(:,fits), lapse(:,fits), Llikelihood(:,fits), LmodelEvidence(:,fits)] = WeightedLinearbiasedLapse_fitter(tsIn(Fit.trialtypes),tpIn(Fit.trialtypes),'InitCond',init,'FitType',Fit.method,[min(tss) max(tss) dt],'N',num2cell(Fit.trialtypes),'LapseSupport',LapseSupport,...
                        'CrossValidation',Fit.CrossValidation,'ModelEvidence',Fit.ModelEvidence,'ObsAct',Fit.ObsAct);
                    
                    % Find likelihood of model on left out condition, if it exists
                    for i = m
                        fitted(i) = any(i == Fit.trialtypes);
                    end
                    if any(~fitted)
                        [notFitLlikelihood(:,fits), notFitLmodelEvidence(:,fits)] = WeigthedLinearLapse_Validator(tsIn(~fitted),tpIn(~fitted),mean(WM(:,fits)),mean(WP(:,fits,1)),mean(B(:,fits)),mean(lapse(:,fits)),...
                            'N',num2cell(m(~fitted)),'LapseSupport',LapseSupport,'ModelEvidence',Fit.ModelEvidence,'FitType',Fit.method,[min(tss) max(tss) dt],'ObsAct',Fit.ObsAct); 
                    else
                        notFitLlikelihood(:,fits) = NaN;
                        notFitLmodelEvidence(:,fits) = NaN;
                    end
                    
                    % Identify lapse trials
                    if fits == Fit.modelUsed;
                        for i = m
                            estimator.type = 'WeightedLinear';
                            estimator.weights = ones(1,Fit.trialtypes(i))/Fit.trialtypes(i);
                            estimator.ObsAct = Fit.ObsAct;
                            estimator.wy = mean(WM(:,fits));
                            [~, ~, loglike, ~, like] = prob_tp_take_ts_wm_wp(tpIn{i}-mean(B(:,fits)),tsIn{i},mean(WM(:,fits)),mean(WP(:,fits,1)),i,'estimator',estimator);
                            loglikeLapse = log(mean(lapse(:,fits))/(LapseSupport(2)-LapseSupport(1)));
                            lapseTrials{i} = log((1-mean(lapse(:,fits)))*like) < loglikeLapse;
                        end
                    end
                    
                    % Recalculate mtp and stdtp
                    for i = m
                        for ii = 1:length(tss)
                            mtp_in(ii,i) = mean(tpIn{i}(tsIn{i} == tss(ii) & ~lapseTrials{i}));
                            stdtp_in(ii,i) = std(tpIn{i}(tsIn{i} == tss(ii) & ~lapseTrials{i}));
                        end
                        errors{i} = tsIn{i}(~lapseTrials{i}) - tpIn{i}(~lapseTrials{i});
                    end
                    
                case 'quad_batch'
                    disp('Fitting measurement averaging model with bias model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estb'
                                init = [0.1 0.06 NaN];      % Default initial search but with estimate of baseline
                                estb(i) = nanmean(tpIn{i}) - nanmean(tsIn{i});
                                init(3) = estb(i);
                            case 'default'
                                init = [0.1 0.06 0];
                        end
                    end
                    dt = Fit.dx;
                    if isfield(varargin{fitnum+1},'batchsize')
                        batchsize = varargin{fitnum+1}.batchsize;
                    else
                        batchsize = 1000000;
                    end
                    LapseSupport = [MinMaxTp(1) max(tss)+MinMaxTp(2)];
                    [WM(:,fits), WP(:,fits,1), B(:,fits), lapse(:,fits), Llikelihood(:,fits), LmodelEvidence(:,fits)] = WeightedLinearbiasedLapse_fitter(tsIn(Fit.trialtypes),tpIn(Fit.trialtypes),'InitCond',init,'FitType',Fit.method,[min(tss) max(tss) dt batchsize],'N',num2cell(Fit.trialtypes),'LapseSupport',LapseSupport,...
                        'CrossValidation',Fit.CrossValidation,'ModelEvidence',Fit.ModelEvidence,'ObsAct',Fit.ObsAct);
                    
                    % Find likelihood of model on left out condition, if it exists
                    for i = m
                        fitted(i) = any(i == Fit.trialtypes);
                    end
                    if any(~fitted)
                        [notFitLlikelihood(:,fits), notFitLmodelEvidence(:,fits)] = WeigthedLinearbiasedLapse_Validator(tsIn(~fitted),tpIn(~fitted),mean(WM(:,fits)),mean(WP(:,fits,1)),mean(B(:,fits)),mean(lapse(:,fits)),...
                            'N',num2cell(m(~fitted)),'LapseSupport',LapseSupport,'ModelEvidence',Fit.ModelEvidence,'FitType',Fit.method,[min(tss) max(tss) dt],'ObsAct',Fit.ObsAct); 
                    else
                        notFitLlikelihood(:,fits) = NaN;
                        notFitLmodelEvidence(:,fits) = NaN;
                    end
                    
                    % Identify lapse trials
                    if fits == Fit.modelUsed
                        for i = m
                            estimator.type = 'WeightedLinear';
                            estimator.weights = ones(1,Fit.trialtypes(i))/Fit.trialtypes(i);
                            estimator.wy = mean(WM(:,fits));
                            [~, ~, loglike, ~, like] = prob_tp_take_ts_wm_wp(tpIn{i}-mean(B(:,fits)),tsIn{i},mean(WM(:,fits)),mean(WP(:,fits,1)),i,'estimator',estimator,'ObsAct',Fit.ObsAct);
                            loglikeLapse = log(mean(lapse(:,fits))/(LapseSupport(2)-LapseSupport(1)));
                            lapseTrials{i} = log((1-mean(lapse(:,fits)))*like) < loglikeLapse;
                        end
                    end
                    
                    % Recalculate mtp and stdtp
                    for i = m
                        for ii = 1:length(tss)
                            mtp_in(ii,i) = mean(tpIn{i}(tsIn{i} == tss(ii) & ~lapseTrials{i}));
                            stdtp_in(ii,i) = std(tpIn{i}(tsIn{i} == tss(ii) & ~lapseTrials{i}));
                        end
                        errors{i} = tsIn{i}(~lapseTrials{i}) - tpIn{i}(~lapseTrials{i});
                    end
                    
                otherwise
                    error('Fitting method not recognized for BLSbiased fitter!')
            end
            G(:,fits) = 1;
            G(:,fits) = 1;
            WM_DRIFT(:,fits) = WM(:,fits);
            W_INT(:,fits) = WM(:,fits)/sqrt(2);
            ALPHA(:,fits) = ones(size(WM,1),1);
            WP(:,fits,2) = WP(:,fits,1);
        
        
        case 'BLS_wm1wm2'
            switch Fit.method
                case 'quad'
                    % Fit the aveMeas model with bias and lapses to the data
                    disp('Fitting measurement BLS_wm1wm2 model with bias and lapses model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estb'
                                init = [0.1 0.06 0.1 0.1/sqrt(2) NaN 0.05];      % Default initial search but with estimate of baseline
                                estb(i) = nanmean(tpIn{i}) - nanmean(tsIn{i});
                                init(5) = estb(i);
                            case 'default'
                                init = [0.1 0.06 0.1 0.1/sqrt(2) 0 0.05];
                        end
                    end
                    dt = Fit.dx;
                    LapseSupport = [MinMaxTp(1) max(tss)+MinMaxTp(2)];
                    [WM(:,fits), WP(:,fits,1), WM_DRIFT(:,fits),...
                     B(:,fits), lapse(:,fits),Llikelihood(:,fits),...
                     LmodelEvidence(:,fits)] = BLS_wm1wm2biasedLapse_fitter(...
                        tsIn(Fit.trialtypes),tpIn(Fit.trialtypes),...
                        'InitCond',init,'FitType',Fit.method,[min(tss) max(tss) dt],...
                        'N',num2cell(Fit.trialtypes),'LapseSupport',LapseSupport,...
                        'CrossValidation',Fit.CrossValidation,...
                        'ModelEvidence',Fit.ModelEvidence,'ObsAct',Fit.ObsAct);
                    
                    % Find likelihood of model on left out condition, if it exists
                    for i = m
                        fitted(i) = any(i == Fit.trialtypes);
                    end
                    if any(~fitted)
                        [notFitLlikelihood(:,fits), notFitLmodelEvidence(:,fits)] = ...
                         BLS_wm1wm2biasedLapse_Validator(...
                            tsIn(~fitted),tpIn(~fitted),...
                            mean(WM(:,fits)),mean(WP(:,fits,1)),...
                            mean(WM_DRIFT(:,fits)),mean(W_INT(:,fits)),...
                            mean(B(:,fits)),mean(lapse(:,fits)),...
                            'N',num2cell(m(~fitted)),...
                            'LapseSupport',LapseSupport,...
                            'ModelEvidence',Fit.ModelEvidence,...
                            'FitType',Fit.method,[min(tss) max(tss) dt],'ObsAct',Fit.ObsAct);
                    else
                        notFitLlikelihood(:,fits) = NaN;
                        notFitLmodelEvidence(:,fits) = NaN;
                    end
                    
                    % Identify lapse trials
                    if fits == Fit.modelUsed
                        for i = m
                            estimator.type = 'BLS_wm1wm2';
                            estimator.wm_drift = mean(WM_DRIFT(:,fits));
                            estimator.ObsAct = Fit.ObsAct;
                            estimator.wy = mean(WP(:,fits,1));
                            [~, ~, loglike, ~, like] = prob_tp_take_ts_wm_wp(...
                                tpIn{i}-mean(B(:,fits)),tsIn{i},...
                                mean(WM(:,fits)),mean(WP(:,fits,1)),i,...
                                'estimator',estimator);
                            loglikeLapse = log(mean(lapse(:,fits))/(LapseSupport(2)-LapseSupport(1)));
                            lapseTrials{i} = log((1-mean(lapse(:,fits)))*like) < loglikeLapse;
                        end
                    end
                    
                    % Recalculate mtp and stdtp
                    for i = m
                        for ii = 1:length(tss)
                            mtp_in(ii,i) = mean(tpIn{i}(tsIn{i} == tss(ii) & ~lapseTrials{i}));
                            stdtp_in(ii,i) = std(tpIn{i}(tsIn{i} == tss(ii) & ~lapseTrials{i}));
                        end
                        errors{i} = tsIn{i}(~lapseTrials{i}) - tpIn{i}(~lapseTrials{i});
                    end
                    
                case 'quad_batch'
                    error('Quad_batch not yet supported!')
                    
                otherwise
                    error('Fitting method not recognized for BLSbiased fitter!')
            end
            G(:,fits) = 1;
            W_INT(:,fits) = nan(size(WM,1),1);
            ALPHA(:,fits) = ones(size(WM,1),1); 
            WP(:,fits,2) = WP(:,fits,1);   
            
        case 'BLS_wm1wm2_wy1wy2'
            switch Fit.method
                case 'quad'
                    % Fit the aveMeas model with bias and lapses to the data
                    disp('Fitting measurement BLS_wm1wm2 model with bias and lapses model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estb'
                                init = [0.1 0.06 0.1 0.1/sqrt(2) NaN 0.05];      % Default initial search but with estimate of baseline
                                estb(i) = nanmean(tpIn{i}) - nanmean(tsIn{i});
                                init(5) = estb(i);
                            case 'default'
                                init = [0.1 0.06 0.1 0.1/sqrt(2) 0 0.05];
                        end
                    end
                    dt = Fit.dx;
                    LapseSupport = [MinMaxTp(1) max(tss)+MinMaxTp(2)];
                    [WM(:,fits), WP(:,fits,1), WP(:,fits,2), WM_DRIFT(:,fits),...
                     B(:,fits), lapse(:,fits),Llikelihood(:,fits),...
                     LmodelEvidence(:,fits)] = BLS_wm1wm2_wy1wy2biasedLapse_fitter(...
                        tsIn(Fit.trialtypes),tpIn(Fit.trialtypes),...
                        'InitCond',init,'FitType',Fit.method,[min(tss) max(tss) dt],...
                        'N',num2cell(Fit.trialtypes),'LapseSupport',LapseSupport,...
                        'CrossValidation',Fit.CrossValidation,...
                        'ModelEvidence',Fit.ModelEvidence,'ObsAct',Fit.ObsAct);
                    
                    % Find likelihood of model on left out condition, if it exists
                    W_INT(:,fits) = nan(size(WM,1),1);
                    for i = m
                        fitted(i) = any(i == Fit.trialtypes);
                    end
                    if any(~fitted)
                        [notFitLlikelihood(:,fits), notFitLmodelEvidence(:,fits)] = ...
                         BLS_wm1wm2_wy1wy2biasedLapse_Validator(...
                            tsIn(~fitted),tpIn(~fitted),...
                            mean(WM(:,fits)),mean(WP(:,fits,1)),...
                            mean(WM_DRIFT(:,fits)),mean(W_INT(:,fits)),...
                            mean(B(:,fits)),mean(lapse(:,fits)),...
                            'N',num2cell(m(~fitted)),...
                            'LapseSupport',LapseSupport,...
                            'ModelEvidence',Fit.ModelEvidence,...
                            'FitType',Fit.method,[min(tss) max(tss) dt],'ObsAct',Fit.ObsAct);
                    else
                        notFitLlikelihood(:,fits) = NaN;
                        notFitLmodelEvidence(:,fits) = NaN;
                    end
                    
                    % Identify lapse trials
                    if fits == Fit.modelUsed
                        for i = m
                            estimator.type = 'BLS_wm1wm2';
                            estimator.wm_drift = mean(WM_DRIFT(:,fits));
                            estimator.w_int = mean(W_INT(:,fits));
                            estimator.ObsAct = Fit.ObsAct;
                            estimator.wy = mean(WP(:,fits,i));
                            [~, ~, loglike, ~, like] = prob_tp_take_ts_wm_wp(...
                                tpIn{i}-mean(B(:,fits)),tsIn{i},...
                                mean(WM(:,fits)),mean(WP(:,fits,i)),i,...
                                'estimator',estimator);
                            loglikeLapse = log(mean(lapse(:,fits))/(LapseSupport(2)-LapseSupport(1)));
                            lapseTrials{i} = log((1-mean(lapse(:,fits)))*like) < loglikeLapse;
                        end
                    end
                    
                    % Recalculate mtp and stdtp
                    for i = m
                        for ii = 1:length(tss)
                            mtp_in(ii,i) = mean(tpIn{i}(tsIn{i} == tss(ii) & ~lapseTrials{i}));
                            stdtp_in(ii,i) = std(tpIn{i}(tsIn{i} == tss(ii) & ~lapseTrials{i}));
                        end
                        errors{i} = tsIn{i}(~lapseTrials{i}) - tpIn{i}(~lapseTrials{i});
                    end
                    
                case 'quad_batch'
                    error('Quad_batch not yet supported!')
                    
                otherwise
                    error('Fitting method not recognized for BLSbiased fitter!')
            end
            G(:,fits) = 1;
            ALPHA(:,fits) = ones(size(WM,1),1); 
            
        case 'EKF'
            switch Fit.method
                case 'quad'
                    % Fit the EKF model with bias and lapses to the data
                    disp('Fitting EKF with bias and lapses model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estb'
                                init = [0.1 0.06 NaN 0.05];      % Default initial search but with estimate of baseline
                                estb(i) = nanmean(tpIn{i}) - nanmean(tsIn{i});
                                init(3) = estb(i);
                            case 'default'
                                init = [0.1 0.06 0 0.05];
                        end
                    end
                    dt = Fit.dx;
                    LapseSupport = [MinMaxTp(1) max(tss)+MinMaxTp(2)];
                    [WM(:,fits), WP(:,fits,1), B(:,fits), lapse(:,fits),...
                        Llikelihood(:,fits), LmodelEvidence(:,fits)] = ...
                        EKF_fitter(...
                        tsIn(Fit.trialtypes),tpIn(Fit.trialtypes),...
                        'InitCond',init,'FitType',Fit.method,[min(tss) max(tss) dt],...
                        'N',num2cell(Fit.trialtypes),'LapseSupport',LapseSupport,...
                        'CrossValidation',Fit.CrossValidation,...
                        'ModelEvidence',Fit.ModelEvidence,'ObsAct',Fit.ObsAct);
                    
                    % Find likelihood of model on left out condition, if it exists
                    for i = m
                        fitted(i) = any(i == Fit.trialtypes);
                    end
                    if any(~fitted)
                        [notFitLlikelihood(:,fits), notFitLmodelEvidence(:,fits)] = EKF_Validator(tsIn(~fitted),tpIn(~fitted),mean(WM(:,fits)),mean(WP(:,fits,1)),mean(B(:,fits)),mean(lapse(:,fits)),...
                            'N',num2cell(m(~fitted)),'LapseSupport',LapseSupport,'ModelEvidence',Fit.ModelEvidence,'FitType',Fit.method,[min(tss) max(tss) dt],'ObsAct',Fit.ObsAct); 
                    else
                        notFitLlikelihood(:,fits) = NaN;
                        notFitLmodelEvidence(:,fits) = NaN;
                    end
                    
                    % Identify lapse trials
                    if fits == Fit.modelUsed
                        estimator.type = 'EKF';
                        estimator.ObsAct = Fit.ObsAct;
                        estimator.wy = mean(WM(:,fits));
                        for i = m
                            [~, ~, loglike, ~, like] = prob_tp_take_ts_wm_wp(tpIn{i}-mean(B(:,fits)),tsIn{i},mean(WM(:,fits)),mean(WP(:,fits,1)),i,'estimator',estimator);
                            loglikeLapse = log(mean(lapse(:,fits))/(LapseSupport(2)-LapseSupport(1)));
                            lapseTrials{i} = log((1-mean(lapse(:,fits)))*like) < loglikeLapse;
                        end
                    end
                    
                    % Recalculate mtp and stdtp
                    for i = m
                        for ii = 1:length(tss)
                            mtp_in(ii,i) = mean(tpIn{i}(tsIn{i} == tss(ii) & ~lapseTrials{i}));
                            stdtp_in(ii,i) = std(tpIn{i}(tsIn{i} == tss(ii) & ~lapseTrials{i}));
                        end
                        errors{i} = tsIn{i}(~lapseTrials{i}) - tpIn{i}(~lapseTrials{i});
                    end
                    
                case 'quad_batch'
                    disp('Fitting BLS with bias model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estb'
                                init = [0.1 0.06 NaN];      % Default initial search but with estimate of baseline
                                estb(i) = nanmean(tpIn{i}) - nanmean(tsIn{i});
                                init(3) = estb(i);
                            case 'default'
                                init = [0.1 0.06 0];
                        end
                    end
                    dt = Fit.dx;
                    if isfield(varargin{fitnum+1},'batchsize')
                        batchsize = varargin{fitnum+1}.batchsize;
                    else
                        batchsize = 1000000;
                    end
                    LapseSupport = [MinMaxTp(1) max(tss)+MinMaxTp(2)];
                    [WM(:,fits), WP(:,fits,1), B(:,fits), lapse(:,fits), Llikelihood(:,fits), LmodelEvidence(:,fits)] = BLSbiasedLapse_fitter(tsIn(Fit.trialtypes),tpIn(Fit.trialtypes),'InitCond',init,'FitType',Fit.method,[min(tss) max(tss) dt batchsize],'N',num2cell(Fit.trialtypes),'LapseSupport',LapseSupport,...
                        'CrossValidation',Fit.CrossValidation,'ModelEvidence',Fit.ModelEvidence,'ObsAct',Fit.ObsAct);
                    
                    % Find likelihood of model on left out condition, if it exists
                    for i = m
                        fitted(i) = any(i == Fit.trialtypes);
                    end
                    if any(~fitted)
                        [notFitLlikelihood(:,fits), notFitLmodelEvidence(:,fits)] = BLSbiasedLapse_Validator(tsIn(~fitted),tpIn(~fitted),mean(WM(:,fits)),mean(WP(:,fits,1)),mean(B(:,fits)),mean(lapse(:,fits)),...
                            'N',num2cell(m(~fitted)),'LapseSupport',LapseSupport,'ModelEvidence',Fit.ModelEvidence,'FitType',Fit.method,[min(tss) max(tss) dt],'ObsAct',Fit.ObsAct); 
                    else
                        notFitLlikelihood(:,fits) = NaN;
                        notFitLmodelEvidence(:,fits) = NaN;
                    end
                    
                    % Identify lapse trials
                    if fits == Fit.modelUsed
                        estimator.type = 'EKF';
                        estimator.ObsAct = Fit.ObsAct;
                        estimator.wy = mean(WM(:,fits));
                        for i = m
                            [~, ~, loglike, ~, like] = prob_tp_take_ts_wm_wp(tpIn{i}-mean(B(:,fits)),tsIn{i},mean(WM(:,fits)),mean(WP(:,fits,1)),i,'estimator',estimator);
                            loglikeLapse = log(mean(lapse(:,fits))/(LapseSupport(2)-LapseSupport(1)));
                            lapseTrials{i} = log((1-mean(lapse(:,fits)))*like) < loglikeLapse;
                        end
                    end
                    
                    % Recalculate mtp and stdtp
                    for i = m
                        for ii = 1:length(tss)
                            mtp_in(ii,i) = mean(tpIn{i}(tsIn{i} == tss(ii) & ~lapseTrials{i}));
                            stdtp_in(ii,i) = std(tpIn{i}(tsIn{i} == tss(ii) & ~lapseTrials{i}));
                        end
                        errors{i} = tsIn{i}(~lapseTrials{i}) - tpIn{i}(~lapseTrials{i});
                    end
                    
                otherwise
                    error('Fitting method not recognized for BLSbiased fitter!')
            end
            G(:,fits) = 1;
            WM_DRIFT(:,fits) = WM(:,fits);
            W_INT(:,fits) = WM(:,fits)/sqrt(2);
            ALPHA(:,fits) = ones(size(WM,1),1);
            WP(:,fits,2) = WP(:,fits,1);
            
        case 'none'
            WM(:,fits) = NaN(m(end),1);
            WP(:,fits,:) = NaN(m(end),1,2);
            B(:,fits) = 0;
            G(:,fits) = 1;
            Llikelihood(:,fits) = NaN(m(end),1);
            estb(:,fits) = NaN(m(end),1);
            lapse(:,fits) = 0;
            for i = 1:length(tsIn)
                lapseTrials{i} = ~ones(size(tsIn{i}));
            end
            G(:,fits) = 1;
            WM_DRIFT(:,fits) = WM(:,fits);
            W_INT(:,fits) = WM(:,fits)/sqrt(2);
            ALPHA(:,fits) = ones(size(WM,1),1);
            
        otherwise
            warning('Fit type argument not recognized, no fit performed.')
            WM(:,fits) = NaN(m(end),1);
            WP(:,fits,:) = NaN(m(end),1,2);
            B(:,fits) = 0;
            G(:,fits) = 1;
            Llikelihood(:,fits) = NaN(m(end),1);
            estb(:,fits) = NaN(m(end),1);
            lapse(:,fits) = 0;
            for i = 1:length(tsIn)
                lapseTrials{i} = ~ones(size(tsIn{i}));
            end
            G(:,fits) = 1;
            WM_DRIFT(:,fits) = WM(:,fits);
            W_INT(:,fits) = WM(:,fits)/sqrt(2);
            ALPHA(:,fits) = ones(size(WM,1),1);
            
    end
end

% Find the parameters to use for the remaining analyses
WMhat = mean(WM(:,Fit.modelUsed));
WPhat = squeeze(mean(WP(:,Fit.modelUsed,:),1));
Bhat = mean(B(:,Fit.modelUsed));
Ghat = mean(G(:,Fit.modelUsed));

% Calculate RMSE, Bias and Variance
for i = m
    rmse(i) = sqrt(mean((errors{i}+Bhat).^2));
    bias(i) = mean((mtp_in(:,i) - Bhat - tss).^2);
    variance(i) = mean(stdtp_in(:,i).^2);
    
    if bootflg
        % Bootstrap bias and variance
        for j = 1:nbootstraps
            es = [];
            for k = 1:length(tss)
                tempP = tpIn{i}(~lapseTrials{i} & tsIn{i} == tss(k));
                inds = ceil(length(tempP)*rand(nsamps,1));
                es = [es; tempP(inds)/Ghat - Bhat - tss(k)];
                tempm(k,:) = mean(tempP(inds));
                tempstd(k,:) = std(tempP(inds));
%                 inds = ceil(length(tp_in{i}{k})*rand(nsamps,1));
%                 tempm(k,:) = mean(tp_in{i}{k}(inds));
%                 tempstd(k,:) = std(tp_in{i}{k}(inds));
            end
            RMSEs(j,i) = sqrt(mean(es.^2));
            BIASs(j,i) = mean((tempm/Ghat - Bhat - tss).^2);
            VARs(j,i) = mean(tempstd.^2);
        end
    else
        BIASs = NaN(1,length(m));
        VARs = NaN(1,length(m));
    end
end

% Generate expected value of aim times based on model fits
switch TAexpectation.method
    case 'none'
        ts_vec = tss;
        ta = nan(length(tss),length(m));
        simbias = nan(1,length(interval_N));
        simv = nan(1,length(interval_N));
    case 'numerical'
        if ~isfield(TAexpectation,'ts_vec') && ~isfield(TAexpectation,'dt')
            ts_vec = tss;
            ts_vec = ts_vec(:);
        elseif ~isfield(TAexpectation,'ts_vec')
            ts_vec = tss(1):TAexpectation.dt:tss(end);
            ts_vec = ts_vec(:);
        else
            ts_vec = TAexpectation.ts_vec;
            ts_vec = ts_vec(:);
        end
        if isfield(TAexpectation,'simtrials')
            simtrials = TAexpectation.simtrials;
        else
            simtrials = 1000;
        end
        for i = m
            switch fittype{Fit.modelUsed}
                case {'BLSbiased','BLS','BLSbiasedLapse'}
                    method_opts.dx = 0.01;
                    if i >= 4
                        ta(:,i) = NaN(size(ts_vec));
                    else
                        ta(:,i) = ta_expectation3(ts_vec,WMhat,i,dt,'Type','BLS','method_options',method_opts,'method','numerical','trials',simtrials,'wp',0,'Support',[min(tss) max(tss)]);
                    end
                    [~, ~, simbias(i), simv(i)] = ta_expectation3(tss',WMhat,i,dt,'method','numerical','trials',simtrials,'wp',WPhat(i),'Support',[min(tss) max(tss)]);
                    
                case 'gBLS'
                    method_opts.dx = 0.01;
                    options.g = Ghat;
                    if i >= 4
                        ta(:,i) = NaN(size(ts_vec));
                    else
                        ta(:,i) = ta_expectation3(ts_vec,WMhat,i,dt,'Type','gBLS','method_options',method_opts,'options',options,'method',TAexpectation.method,'Support',[min(tss) max(tss)]);
                    end
                    
                    if i == 1
                        [~, ~, simbias(1), simv(1)] = ta_expectation3(tss',WMhat,1,dt,'method','numerical','trials',1000,'wp',WPhat(i),'Support',[min(tss) max(tss)]);
                    elseif i == 2
                        [~, ~, simbias(2), simv(2)] = ta_expectation3(tss',WMhat,2,dt,'method','numerical','trials',1000,'wp',WPhat(i),'Support',[min(tss) max(tss)]);
                    elseif i > 2
                        [~, ~, simbias(i), simv(i)] = ta_expectation3(tss',WMhat,i,dt,'method','numerical','trials',1000,'wp',WPhat(i),'Support',[min(tss) max(tss)]);
                    end
                case {'ObsAct','ObsActBiased','ObsActBiasedLapse'}
                    method_opts.dx = 0.01;
                    if i >= 4
                        [ta(:,i), ta_std(:,i)] = NaN(size(ts_vec));
                    else
                        [ta(:,i), ta_std(:,i)] = ta_expectation3(ts_vec,WMhat,i,dt,'Type','ObsAct','method_options',method_opts,'method','numerical','trials',simtrials,'wp',WPhat(i),'Support',[min(tss) max(tss)]);
                    end
                    [~, ~, simbias(i), simv(i)] = ta_expectation3(tss',WMhat,i,dt,'Type','ObsAct','method','numerical','trials',simtrials,'wp',WPhat(i),'Support',[min(tss) max(tss)]);
                
                case 'aveMeasurements'
                    method_opts.dx = 0.01;
                    ta(:,i) = ta_expectation3(ts_vec,WMhat,i,dt,'Type','aveMeasurements','method_options',method_opts,'method','numerical','trials',simtrials,'wp',0,'Support',[min(tss) max(tss)]);
                    [~, ~, simbias(i), simv(i), simrmse(i)] = ta_expectation3(tss',WMhat,i,dt,'Type','aveMeasurements','method_options',method_opts,'method','numerical','trials',simtrials,'wp',WPhat(i),'Support',[min(tss) max(tss)]);
                    
                case 'MAP'
                    method_opts.dx = 0.01;
                    ta(:,i) = ta_expectation3(ts_vec,WMhat,i,dt,'Type','MAP','method_options',method_opts,'method','numerical','trials',simtrials,'wp',0,'Support',[min(tss) max(tss)]);
                    [~, ~, simbias(i), simv(i), simrmse(i)] = ta_expectation3(tss',WMhat,i,dt,'Type','MAP','method_options',method_opts,'method','numerical','trials',simtrials,'wp',WPhat(i),'Support',[min(tss) max(tss)]);
                    
                case 'WeightedLinear'
                    method_opts.dx = 0.01;
                    ta(:,i) = ta_expectation3(ts_vec,WMhat,i,dt,'Type','WeightedLinear','method_options',method_opts,'method','numerical','trials',simtrials,'wp',0,'Support',[min(tss) max(tss)]);
                    [~, ~, simbias(i), simv(i), simrmse(i)] = ta_expectation3(tss',WMhat,i,dt,'Type','WeightedLinear','method_options',method_opts,'method','numerical','trials',simtrials,'wp',WPhat(i),'Support',[min(tss) max(tss)]);
                    
                    
                case 'SubOptMemBias'
                    method_opts.dx = 0.01;
                    estimator.type = 'SubOptMemBias';
                    WMhat = mean(WM(:,Fit.modelUsed));
                    estimator.wm_drift = mean(WM_DRIFT(:,Fit.modelUsed));
                    estimator.w_int = mean(W_INT(:,Fit.modelUsed));
                    estimator.ObsAct = Fit.ObsAct;
                    estimator.wy = mean(WP(:,Fit.modelUsed,i));
                    ta(:,i) = ta_expectation3(ts_vec,WMhat,i,dt,'Type','N/A','method_options',method_opts,'method','numerical','trials',simtrials,'wp',0,'Support',[min(tss) max(tss)],'estimator',estimator);
                    [~, ~, simbias(i), simv(i), simrmse(i)] = ta_expectation3(tss',WMhat,i,dt,'Type','N/A','method_options',method_opts,'method','numerical','trials',simtrials,'wp',WPhat(i),'Support',[min(tss) max(tss)],'estimator',estimator);
                    
                    
                case 'NestedModel'
                    method_opts.dx = 0.01;
                    estimator.type = 'NestedModel';
                    WMhat = mean(WM(:,Fit.modelUsed));
                    estimator.wm_drift = mean(WM_DRIFT(:,Fit.modelUsed));
                    estimator.w_int = mean(W_INT(:,Fit.modelUsed));
                    estimator.alpha = mean(ALPHA(:,Fit.modelUsed));
                    estimator.ObsAct = Fit.ObsAct;
                    estimator.wy = mean(WP(:,Fit.modelUsed,i));
                    ta(:,i) = ta_expectation3(ts_vec,WMhat,i,dt,'Type','N/A','method_options',method_opts,'method','numerical','trials',simtrials,'wp',0,'Support',[min(tss) max(tss)],'estimator',estimator);
                    [~, ~, simbias(i), simv(i), simrmse(i)] = ta_expectation3(tss',WMhat,i,dt,'Type','N/A','method_options',method_opts,'method','numerical','trials',simtrials,'wp',WPhat(i),'Support',[min(tss) max(tss)],'estimator',estimator);
                    
                case {'BLS_wm1wm2','BLS_wm1wm2_wy1wy2'}
                    method_opts.dx = 0.01;
                    estimator.type = 'BLS_wm1wm2';
                    WMhat = mean(WM(:,Fit.modelUsed));
                    estimator.wm_drift = mean(WM_DRIFT(:,Fit.modelUsed));
                    estimator.ObsAct = Fit.ObsAct;
                    estimator.wy = mean(WP(:,Fit.modelUsed,i));
                    ta(:,i) = ta_expectation3(ts_vec,WMhat,i,dt,'Type','N/A','method_options',method_opts,'method','numerical','trials',simtrials,'wp',0,'Support',[min(tss) max(tss)],'estimator',estimator);
                    [~, ~, simbias(i), simv(i), simrmse(i)] = ta_expectation3(tss',WMhat,i,dt,'Type','N/A','method_options',method_opts,'method','numerical','trials',simtrials,'wp',WPhat(i),'Support',[min(tss) max(tss)],'estimator',estimator);
                    
                case {'EKF'}
                    method_opts.dx = 0.01;
                    estimator.type = 'EKF';
                    WMhat = mean(WM(:,Fit.modelUsed));
                    estimator.wm_drift = mean(WM_DRIFT(:,Fit.modelUsed));
                    estimator.ObsAct = Fit.ObsAct;
                    estimator.wy = mean(WP(:,Fit.modelUsed,i));
                    ta(:,i) = ta_expectation3(ts_vec,WMhat,i,dt,'Type','N/A','method_options',method_opts,'method','numerical','trials',simtrials,'wp',0,'Support',[min(tss) max(tss)],'estimator',estimator);
                    [~, ~, simbias(i), simv(i), simrmse(i)] = ta_expectation3(tss',WMhat,i,dt,'Type','N/A','method_options',method_opts,'method','numerical','trials',simtrials,'wp',WPhat(i),'Support',[min(tss) max(tss)],'estimator',estimator);
                    
                otherwise
                    ta = NaN;
                    simbias(i) = NaN;
                    simv(i) = NaN;
                    simrmse(i) = NaN;
            end
        end
        
    case 'analytical'
        if ~isfield(TAexpectation,'ts_vec') && ~isfield(TAexpectation,'dt')
            ts_vec = tss;
            ts_vec = ts_vec(:);
        elseif ~isfield(TAexpectation,'ts_vec')
            ts_vec = tss(1):TAexpectation.dt:tss(end);
            ts_vec = ts_vec(:);
        else
            ts_vec = TAexpectation.ts_vec;
            ts_vec = ts_vec(:);
        end
        if isfield(TAexpectation,'simtrials')
            simtrials = TAexpectation.simtrials;
        else
            simtrials = 1000;
        end
        for i = m
            switch fittype{Fit.modelUsed}
                case {'BLSbiased','BLS'}
                    method_opts.dx = 0.01;
                    if i >= 4
                        ta(:,i) = NaN(size(ts_vec));
                    else
                        ta(:,i) = ta_expectation3(ts_vec,WMhat,i,dt,'Type','BLS','method_options',method_opts,'method','analytical','Support',[min(tss) max(tss)]);
                    end
                    [~, ~, simbias(i), simv(i)] = ta_expectation3(tss,WMhat,i,dt,'method','numerical','trials',simtrials,'wp',WPhat(i),'Support',[min(tss) max(tss)]);
                    
                case 'gBLS'
                    method_opts.dx = 0.01;
                    options.g = Ghat;
                    if i >= 4
                        ta(:,i) = NaN(size(ts_vec));
                    else
                        ta(:,i) = ta_expectation3(ts_vec,WMhat,i,dt,'Type','gBLS','method_options',method_opts,'options',options,'method',TAexpectation.method,'Support',[min(tss) max(tss)]);
                    end
                    
                    if i == 1
                        [~, ~, simbias(1), simv(1)] = ta_expectation3(tss',WMhat,1,dt,'method','numerical','trials',1000,'wp',WPhat(i),'Support',[min(tss) max(tss)]);
                    elseif i == 2
                        [~, ~, simbias(2), simv(2)] = ta_expectation3(tss',WMhat,2,dt,'method','numerical','trials',1000,'wp',WPhat(i),'Support',[min(tss) max(tss)]);
                    elseif i > 2
                        [~, ~, simbias(i), simv(i)] = ta_expectation3(tss',WMhat,i,dt,'method','numerical','trials',1000,'wp',WPhat(i),'Support',[min(tss) max(tss)]);
                    end
                case 'ObsAct'
                    error('Not yet supported')
                    ta(:,i) = ta_expectation3(ts_vec,WMhat,1,'Type','ObsAct','wp',wp(i),'Support',[min(tss) max(tss)]);
                    
                    if i == 2
                        ta2 = ta_expectation(ts_vec,wm(1),2,'Type','ObsAct','wp',wp(1),'Support',[min(tss) max(tss)]);
                    elseif i > 2
                        error('Not yet supported')
                    end
                    
                case 'aveMeasurements'
                    method_opts.dx = 0.01;
                    if i >= 4
                        ta(:,i) = NaN(size(ts_vec));
                    else
                        ta(:,i) = ta_expectation3(ts_vec,WMhat,i,dt,'Type','aveMeasurements','method_options',method_opts,'method','analytical','Support',[min(tss) max(tss)]);
                    end
                    [~, ~, simbias(i), simv(i), simrmse(i)] = ta_expectation3(tss,WMhat,i,dt,'Type','aveMeasurements','method','numerical','trials',simtrials,'wp',WPhat(i),'Support',[min(tss) max(tss)]);
                    
                otherwise
                    ta = NaN;
                    simbias(i) = NaN;
                    simv(i) = NaN;
                    
            end
        end
        
    otherwise
        error(['Model expectation method ' TAexpectation.method ' not recognized!']);
end

%% PLOTTING
switch Plot
    case {'All','all'}
        titles = PlotOpts.titles;
        RelativeFigSize = PlotOpts.RelativeFigSize;
        colors = PlotOpts.colors;
        
        % Dependence on sample time
        scrsz = get(groot,'ScreenSize');
        figure('Name',[d.sname ' dependence on sample time'],'Position',[scrsz(3) scrsz(4) scrsz(3) scrsz(4)].*RelativeFigSize)
        maxrmse = max(max(rmse));
        plotind = 1;
        allts = [];
        alltp = [];
        for i = m
            for ii = 1:length(ts_in{i})
                allts = [allts; ts_in{i}{ii}];
                alltp = [alltp; tp_in{i}{ii}];
            end
        end
        ax = [min(allts)-100 max(allts)+100 min(alltp)-100 max(alltp)+100];
        for i = m
            subplot(1,length(m),plotind)
            plotind = plotind+1;
            plot(tsIn{i},tpIn{i},'o','Color',colors(i,:)+(1 - colors(i,:))/1.5)
            hold all
            plot(tsIn{i}(lapseTrials{i}),tpIn{i}(lapseTrials{i}),'.','Color',[0 0 0])
%             for ii = 1:length(ts_in{i})
%                 plot(ts_in{i}{ii},tp_in{i}{ii},'.','Color',colors(i,:)+(1 - colors(i,:))/1.5)
%                 hold all
%             end
            plot(tss-200:tss(end)+200,tss(1)-200:tss(end)+200,'k')
            text(tss(end),tss(1)-100,['p = ' num2str(pval(i))]);
            axis(ax)
            xlabel('Sample time (ms)')
            ylabel('Production time (ms)')
            title(titles{i});
        end
        
        plotind = 1;
        for i = m
            subplot(1,length(m),plotind)
            plotind = plotind+1;
            errorbar(tss,mtp_in(:,i),stdtp_in(:,i),'.','Color',colors(i,:),'MarkerSize',20,'LineWidth',2)
            hold all
            if ~strcmp('none',fittype{Fit.modelUsed}) && any(i == interval_N)
                plot(ts_vec,ta(:,i)+Bhat,'Color',colors(i,:))
            end
            axis(ax)
            plot(tss(1)-200:tss(end)+200,tss(1)-200:tss(end)+200,'k')
        end
        xlabel('Sample time (ms)')
        ylabel('Mean production time (ms)')
        
        figure('Name',[d.sname ' bias vs. sqrt(variance)'])
        for i = m
            h(i) = plot(0:0.01:rmse(i),sqrt(rmse(i)^2-(0:0.01:rmse(i)).^2),'Color',colors(i,:));
            hold on
            if bootflg
                plot(sqrt(BIASs(:,i)),sqrt(VARs(:,i)),'.','Color',colors(i,:) + 0.7*(colors(i,:)==0))
            end
            if ~strcmp('none',fittype{Fit.modelUsed})
                h(i) = plot(sqrt(simbias(i)),sqrt(simv(i)),'o','Color',colors(i,:));
            end
            h(i) = plot(sqrt(bias(i)),sqrt(variance(i)),'.','Color',colors(i,:),'MarkerSize',20);
        end
        axis([0 maxrmse+50 0 maxrmse+50])
        axis square
        xlabel('bias')
        ylabel('sqrt(variance)')
        legend(h,titles)
        
    case {'No','no','NO','N','n','none','None',0}
    otherwise
        error('Plot option not recognized!');
        
end
        
%% Saving
switch Save
    case {'default','YES','Yes','yes','y','Y',1}
        % Save the results to the default location
        filename = [d.spath '/' d.sname '_' d.projname '_PoolAnalysisResults' datestr(now,'yyyymmdd')];
        disp(['Saving results to ' filename])
        save(filename);
        
    case {'No','no','NO','n','N',0}
        % Don't save the results
    otherwise
        save(Save);     % Save to file specified by the Save string
        
end
