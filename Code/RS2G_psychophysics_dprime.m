function [dprime, mtp_in, stdtp_in, pval, weber, tsIn, tpIn, trialsIn, ts_in, tp_in, Trials_sorted] = RS2G_psychophysics_dprime(d,varargin)
%% RS4GPoolAnalysis
%
%   [mtp_in, stdtp_in, bias, variance, rmse, wm, wp, B, pval, weber, tsIn, tpIn, trialsIn, ts_in, tp_in, Trials_sorted, estb, llikelihood, ta, ta2] = RS4GPoolAnalysis(varargin)   
%
%%

%% Defaults
PlotOpts_default.titles = {'RS1G','RS2G'};
PlotOpts_default.RelativeFigSize = [1/5 1/2 3/5 3/5];
PlotOpts_default.colors = [0 0 1; 1 0 0; 0.6 0.6 0.6; 0 0 0];
PlotOpts_default.nbins = 20;

%% Parse input
Parser = inputParser;

addRequired(Parser,'d')     % Data structure
addParameter(Parser,'runs',NaN)     % Runs to analyze
addParameter(Parser,'interval_N',1:2)    % Trial types to analyze
addParameter(Parser,'sampleIntervals',600:100:1000)   % Sample intervals to measure sensitivity to
addParameter(Parser,'outlier',Inf)      % Number of standard deviations away from the mean to exclude data as outlier
addParameter(Parser,'trialWin',[1 Inf]) % Trials in each run to analyze ([start# end#])
addParameter(Parser,'Plot','none')      % Data to plot
addParameter(Parser,'PlotOpts',PlotOpts_default)    % Plotting options
addParameter(Parser,'Save','No')           % Saving options
addParameter(Parser,'MinMaxTp',[-Inf Inf])      % Minimum and maximum values of t_p to keep
addParameter(Parser,'ConflictType','equal')     % For experiments with cue conflict
addParameter(Parser,'DiffTolerance',2/60)       % Tolerance for difference in sample times before calling it conflict

parse(Parser,d,varargin{:})

d = Parser.Results.d;
runs = Parser.Results.runs;
interval_N = Parser.Results.interval_N;
sampleIntervals = Parser.Results.sampleIntervals;
outlier = Parser.Results.outlier;
trialWin = Parser.Results.trialWin;
Plot = Parser.Results.Plot;
PlotOpts = Parser.Results.PlotOpts;
Save = Parser.Results.Save;
MinMaxTp = Parser.Results.MinMaxTp;
ConflictType = Parser.Results.ConflictType;
DiffTolerance = Parser.Results.DiffTolerance;

% Check to see if run information was provided
if isnan(runs)
    runs = 1:d.runs;            % Defaults to pool data from all runs
end

% Set m to the length of interval_N
m = 1:length(interval_N);

% Set nbins
nbins = PlotOpts.nbins;

%% Analyze data

for i = m
    % Grab the appropriate data
    [ts1{i}, ts2{i}, ts{i}, tp{i}, ~, ~, Trials{i}, correct{i}, ~, ~, ~, N{i}] = RS2G_psychophysics_pooldata(d,'runs',runs,'interval_N',interval_N(i),'trialWin',trialWin);
    
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
    tss{i} = unique(ts{i});
    
    % Estimate Weber fraction
    disp('Estimating Weber fraction')
    if isempty(tp{i})
        weber(i) = 0.1;
    elseif sum(tp{i} > 0)
        fun = @(w)(sum(w.^2*tss{i}.^2) - var(tp{i}(tp{i} > 0 & tp{i} >= MinMaxTp(1) & tp{i} <= tss{i}(end)+MinMaxTp(2))));
        weber(i) = lsqnonlin(fun,0.1);
    else
        weber(i) = NaN;
    end
    
    % Sort data by sample time
    errors{i} = [];
    tsIn{i} = [];
    tpIn{i} = [];
    trialsIn{i} = [];
    for ii = 1:length(tss{i})
        ts_sorted{i}{ii} = ts{i}(ts{i} == tss{i}(ii));
        tp_sorted{i}{ii} = tp{i}(ts{i} == tss{i}(ii));
        Trials_sorted{i}{ii} = Trials{i}(ts{i} == tss{i}(ii));
        
        mtp(ii,i) = nanmean(tp_sorted{i}{ii}(tp_sorted{i}{ii} >= MinMaxTp(1) & tp_sorted{i}{ii} <= tss{i}(ii)+MinMaxTp(2)));
        
        stdtp(ii,i) = weber(i)*tss{i}(ii);
        
        % Find outliers and recalculate mean and variance
        ins{i}{ii} = find(abs(tp_sorted{i}{ii} - mtp(ii,i)) < outlier*stdtp(ii,i) & tp_sorted{i}{ii} >= MinMaxTp(1) & tp_sorted{i}{ii} <= tss{i}(ii)+MinMaxTp(2));
        
        ts_in{i}{ii} = ts_sorted{i}{ii}(ins{i}{ii});
        tp_in{i}{ii} = tp_sorted{i}{ii}(ins{i}{ii});
        trials_in{i}{ii} = Trials_sorted{i}{ii}(ins{i}{ii});
        tsIn{i} = [tsIn{i}; ts_sorted{i}{ii}(ins{i}{ii})];
        tpIn{i} = [tpIn{i}; tp_sorted{i}{ii}(ins{i}{ii})];
        trialsIn{i} = [trialsIn{i}; Trials_sorted{i}{ii}(ins{i}{ii})];
        
        mtp_in(ii,i) = nanmean(tp_sorted{i}{ii}(ins{i}{ii}));
        
        stdtp_in(ii,i) = nanstd(tp_sorted{i}{ii}(ins{i}{ii}),1);
        
        errors{i} = [errors{i}; tss{i}(ii) - tp_in{i}{ii}];
    end
    
    % Find the shortest and longest production intervals
    tpMinMax = [min(tpIn{i}) max(tpIn{i})];
    
    % Test if production times for longest interval are significantly longer than shortest interval
    tpshort = tp_in{i}{1};
    tplong = tp_in{i}{end};
    
    [h, pval(i)] = ttest2(tplong,tpshort,'tail','right');
    
    % Measure d' for each combination of sample intervals
    for j = 1:length(sampleIntervals)
        sampIntInd(j) = find(tss{i} == sampleIntervals(j));
    end
    [MUa, MUb] = meshgrid(mtp_in(sampIntInd,i));
    [SIGa, SIGb] = meshgrid(stdtp_in(sampIntInd,i));
    dprime(:,:,i) = ( MUa-MUb ) ./ sqrt( (SIGa + SIGb)/2 );
    
end         % interval_N loop



%% PLOTTING
switch Plot
    case {'All','all'}
        titles = PlotOpts.titles;
        RelativeFigSize = PlotOpts.RelativeFigSize;
        colors = PlotOpts.colors;
        
        % Distribution of responses
        scrsz = get(groot,'ScreenSize');
        bins = linspace(tpMinMax(1),tpMinMax(2),nbins);
        for i = m
            figure('Name',[d.sname ' interval_N = ' num2str(interval_N(i))],'Position',[scrsz(3) scrsz(4) scrsz(3) scrsz(4)].*RelativeFigSize)
            plotind = 0;
            for j = 1:length(sampleIntervals)
                Nj = hist(tp_in{i}{sampIntInd(j)},bins);
                Pj = Nj/sum(Nj);
                for k = 1:length(sampleIntervals)
                    plotind = plotind+1;
                    if j <= k
                        Nk = hist(tp_in{i}{sampIntInd(k)},bins);
                        Pk = Nk/sum(Nk);
                        subplot(length(sampleIntervals),length(sampleIntervals),plotind)
                        bar(bins,Pj,'FaceColor',colors(1,:),'EdgeColor',colors(1,:));
                        hold on
                        bar(bins,Pk,'FaceColor',colors(2,:),'EdgeColor',colors(2,:));
                        plot(bins,normpdf(bins,mtp_in(sampIntInd(j),i),stdtp_in(sampIntInd(j),i))*(bins(2)-bins(1)),'Color',colors(1,:),'LineWidth',2)
                        plot(bins,normpdf(bins,mtp_in(sampIntInd(k),i),stdtp_in(sampIntInd(k),i))*(bins(2)-bins(1)),'Color',colors(2,:),'LineWidth',2)
                        maxP = max([Pj(:); Pk(:)]);
                        plot(mtp_in(sampIntInd(j),i)*ones(1,2),[0 1.2*maxP],'k','LineWidth',2)
                        plot(mtp_in(sampIntInd(k),i)*ones(1,2),[0 1.2*maxP],'k','LineWidth',2)
                        plot([mtp_in(sampIntInd(j),i) mtp_in(sampIntInd(k),i)],1.2*maxP*ones(1,2),'k','LineWidth',2)
                        text((mtp_in(sampIntInd(j),i)+mtp_in(sampIntInd(k),i))/2,1.25*maxP,['d'' = ' num2str(dprime(j,k,i))],'Interpreter','latex')
                        xlabel('t_p (ms)')
                        ylabel('p(t_p)')
                    end
                end
            end    
        end
        
        figure('Name',[d.sname ' d'' matrix'])
        for i = m
            subplot(1,length(interval_N),i)
            imagesc(sampleIntervals,sampleIntervals,dprime(:,:,i))
            colormap gray
            xlabel('t_s (ms)')
            ylabel('t_s (ms)')
            axis square
            title(['interval_N = ' num2str(interval_N(i))])
            colorbar
        end
        
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
