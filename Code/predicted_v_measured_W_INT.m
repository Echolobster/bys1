function [BIAS, sqrtVAR, SimBiasBLS, SimVarBLS, SimBiasAve, SimVarAve,...
    deltaBV, deltaBVBLS, deltaBVAve] = predicted_v_measured_W_INT(...
    slist,varargin)
%% predicted_v_actualBiasVariance
%
%   [biases, variances, biasBLS, varBLS, biasAve, varAve, deltaBV,
%   deltaBVBLS, deltaBVAve] = predicted_v_actualBiasVariance(list)
%
%       Computes the predicted bias and variance given fits to the BLS and
%       averaging models and compares agains the observed bias and
%       variance.
%
%%

%% Defaults
PlotOpts_default.colors = [0 0 1; 1 0 0];

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'slist')     % List of subjects to analyze
addParameter(Parser,'N',2)      % Maximum number of sets
addParameter(Parser,'tss',600:100:1000)     % sample times for experiment
addParameter(Parser,'simulationN',100000)    % Number of trials per simulation
addParameter(Parser,'CommonFileName','_BLSbiasedFitResults20170412')
addParameter(Parser,'vecN',100)
addParameter(Parser,'Plot','Yes')
addParameter(Parser,'PlotOpts',PlotOpts_default)

parse(Parser,slist,varargin{:})

slist = Parser.Results.slist;
N = Parser.Results.N;
tss = Parser.Results.tss;
simulationN = Parser.Results.simulationN;
CommonFileName = Parser.Results.CommonFileName;
vecN = Parser.Results.vecN;
Plot = Parser.Results.Plot;
PlotOpts = Parser.Results.PlotOpts;

if length(simulationN) == 1
    simulationN = simulationN*ones(length(slist),1);
end

%% Load model fits and observed bias and variance for each subject

for i = 1:length(slist)
    load([slist{i} CommonFileName],'WM','WP','B','lapse','bias','BIASs',...
        'variance','VARs','rmse','RMSEs','Fit','Llikelihood','mtp_in',...
        'notFitLlikelihood','stdtp_in','WM_DRIFT','W_INT','ALPHA',...
        'tsIn','lapseTrials')
    
    % BIAS/VAR
    BIAS(i,:) = sqrt(bias);
    sqrtVAR(i,:) = sqrt(variance);
    RMSE(i,:) = rmse;
    
    % Parameters fit to data
    wm(i,:) = mean(WM,1);
    wp(i,:) = mean(WP(:,:,1),1);
    wp2(i,:) = mean(WP(:,:,2),1);
    b(i,:) = mean(B,1);
    Lapse(i,:) = mean(lapse,1);
    if exist('WM_DRIFT','var')
        if size(WM_DRIFT,2) == 1
            wm_drift(i,:) = [mean(WM_DRIFT,1) NaN];
        else
            wm_drift(i,:) = mean(WM_DRIFT,1);
        end
    else
        wm_drift(i,:) = [NaN NaN];
    end
    if exist('W_INT','var')
        if size(W_INT,2) == 1
            w_int(i,:) = [mean(W_INT,1) NaN];
        else
            w_int(i,:) = mean(W_INT,1);
        end
    else
        w_int(i,:) = [NaN NaN];
    end
    if exist('ALPHA','var')
        if size(ALPHA,2) == 1
            alpha(i,:) = [mean(ALPHA,1) NaN];
        else
            alpha(i,:) = mean(ALPHA,1);
        end
    else
        alpha(i,:) = [NaN NaN];
    end
    
    % Bootleg BIAS/VAR
    Bs{i} = sqrt(BIASs);
    Vs{i} = sqrt(VARs);
%     RMSError{i} = RMSEs;
    RMSError{i} = sqrt(BIASs+VARs);
    
    % Model info
    ModelUsed = Fit.modelUsed;
    Models{i} = Fit.fittype;
%     LLmodels{i} = -notFitLlikelihood;
    LLmodels{i} = -Llikelihood(1:end-1,:);  % TODO fix Llikelihood to be per trial
    
    % Mean and standard deviation of production times
    MTP(:,:,i) = mtp_in;
    STDTP(:,:,i) = stdtp_in;
    
    
    % Calculate simulationN
    if any(isnan(simulationN))
        simulationN(i) = floor(...
            (numel(tsIn{1}(~lapseTrials{1})) + numel(tsIn{2}(~lapseTrials{2})))/...
            2/length(tss));      % average trials/ts
    end
end

wm_mem = wm_drift;

if length(unique(ModelUsed)) > 1
    error('Subjects were fit to different models!')
else
    ModelUsed = unique(ModelUsed);
end



%% Determine the best w_int for each subject's wm and wm_mem
wVec = nan(vecN,length(slist));
RMSEopt = nan(length(slist),1);
RMSEtest = nan(vecN,length(slist));
w_intOPT = nan(length(slist),1);


for i = 1:length(slist)
    
    % Parameters
    estimator.type = Models{i}{ModelUsed};
    estimator.wm_drift = wm(i,ModelUsed);%wm_mem(i,ModelUsed);
    estimator.w_int = w_int(i,ModelUsed);
    wVec(:,i) = linspace(wm(i,ModelUsed)/sqrt(2)/4,...
        1.1*max([wm(i,ModelUsed) wm_mem(i,ModelUsed)]),vecN);
    
    % Find w_int
    [w_intOPT(i), RMSEopt(i), RMSEtest(:,i)] = optimalW_INT(wm(i,ModelUsed),...
        wm_mem(i,ModelUsed),tss,'estimator',estimator,...
        'wVec',wVec(:,i),'trials',simulationN(i),'wp',wp(i,ModelUsed));
    
end

%% Plotting
switch Plot
    case {'Yes','yes','y','Y','YES'}
        
        %% Set up colors
        colors =[ 0.3010    0.7450    0.9330;...
            0.8500    0.3250    0.0980;...
            0.9290    0.6940    0.1250;...
            0.4940    0.1840    0.5560;...
            0.4660    0.6740    0.1880;...
            0         0.4470    0.7410;...
            0.6350    0.0780    0.1840;...
            0         0         1.0000;...
            0         0.5000    0;...
            1.0000    0         0];
        set(groot,'defaultAxesColorOrder',colors)
        
        
        %% w_int optimal vs w_int observed
        figure('Name','w_{int}^* vs w_{int}')
        for i = 1:length(slist)
            h(i) = plot(w_intOPT(i),w_int(i,ModelUsed),'o','Color',colors(i,:),...
                'MarkerFaceColor',colors(i,:));
            hold on
        end
        axis([0 max([w_intOPT(:); w_int(:,ModelUsed(1)); wm(:,ModelUsed(1))]) ...
            0 max([w_intOPT(:); w_int(:,ModelUsed(1)); wm(:,ModelUsed(1))])])
        axis square
        plotUnity;
        xlabel('w_{int}^*')
        ylabel('w_{int}')
        legend(h,slist)
        mymakeaxis(gca);
        
        %% RMSE as a function of w_int
        figure('Name','RMSE as a function of w_{int}')
        for i = 1:length(slist)
            plot(wVec(:,i),RMSEtest(:,i),'Color',colors(i,:))
            hold on
            lineProps.LineStyle = '--';
            lineProps.Color = colors(i,:);
            plotVertical(w_int(i,ModelUsed),...
                'MinMax',[0 RMSE(i,2)],...
                'lineProperties',lineProps);
            lineProps.LineStyle = '-';
            lineProps.Color = colors(i,:);
            plotVertical(w_intOPT(i),...
                'MinMax',[0 RMSEopt(i)],...
                'lineProperties',lineProps);
        end
        axis tight
        xlabel('w_{int}')
        ylabel('RMSE')
        mymakeaxis(gca)
        
end
