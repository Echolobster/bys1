function [BIAS, sqrtVAR, SimBiasBLS, SimVarBLS, SimBiasAve, SimVarAve, deltaBV, deltaBVBLS, deltaBVAve] = predicted_v_actualBiasVariance(slist,varargin)
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

% Defaults
TAexpectation.methodopts.dx = 0.01;
TAexpectation.dt = 10;
PlotOpts_default.colors = [1 0 0; 0 0 1];

% Parse inputs
Parser = inputParser;

addRequired(Parser,'slist')     % List of subjects to analyze
addParameter(Parser,'N',2)      % Maximum number of sets
addParameter(Parser,'tss',600:100:1000)     % sample times for experiment
addParameter(Parser,'simulationN',10000)    % Number of trials per simulation
addParameter(Parser,'CommonFileName','_BLSbiasedFitResults20150809')
addParameter(Parser,'TAexpectation',TAexpectation_default)  % For controlling the calculation of the expected value of aim times under a model
addParameter(Parser,'Plot','Yes')
addParameter(Parser,'PlotOpts',PlotOpts_default)

parse(Parser,slist,varargin{:})

slist = Parser.Results.slist;
N = Parser.Results.N;
tss = Parser.Results.tss;
simulationN = Parser.Results.simulationN;
CommonFileName = Parser.Results.CommonFileName;
TAexpectation = Parser.Results.TAexpectation;
Plot = Parser.Results.Plot;
PlotOpts = Parser.Results.PlotOpts;

%% Load model fits and observed bias and variance for each subject

for i = 1:length(slist)
    load([slist{i} CommonFileName],'WM','WP','bias','variance')
    
    BIAS(i,:) = sqrt(bias);
    sqrtVAR(i,:) = sqrt(variance);
    wm(i) = WM;
    wp(i) = WP;
    
end

%% Simulate the BLS and averaging models for the parameters fit to each subject

for i = 1:length(slist)
    for j = 1:N
        [~, ~, simbias, simv] = ta_expectation3(tss',wm(i),j,TAexpectation.dt,'method','numerical','trials',simulationN,'wp',WP,'Support',[min(tss) max(tss)]);
        SimBiasBLS(i,j) = sqrt(simbias);
        SimVarBLS(i,j) = sqrt(simv);
        
        [~, ~, simbias, simv] = ta_expectation3(tss',wm(i),j,TAexpectation.dt,'method','numerical','trials',simulationN,'wp',WP,'Support',[min(tss) max(tss)],'Type','aveMeasurments');
        SimBiasAve(i,j) = sqrt(simbias);
        SimVarAve(i,j) = sqrt(simv);
    end
end

%% Find the differences between the bias and variances of data and models
deltaBV(:,1) = diff(BIAS,[],2);
deltaBV(:,2) = diff(sqrtVAR,[],2);

deltaBVBLS(:,1) = diff(SimBiasBLS,[],2);
deltaBV(:,2) = diff(SimVarBLS,[],2);

deltaBVAve(:,1) = diff(SimBiasAve,[],2);
deltaBVAve(:,2) = diff(SimVarAve,[],2);

%% Find the dot product between model predicted changes and the actual changes

for i = 1:length(slist)
    dp(i,1) = deltaBV'*deltaBVBLS;
    dp(i,2) = deltaBV'*deltaBVAve;
end


%% Ploting
switch Plot
    case {'Yes','yes','y','Y','YES'}
        figure('Name','Actual vs BLS prediction')
        subplot(1,2,1)
        for n = 1:N
            plot(BIAS(:,n),SimBiasBLS(:,n),'.','Color',PlotOpts.colors(n,:))
            hold on
        end
        xlabel('Bias (ms)')
        ylabel('Expected bias (ms)')
        
        subplot(1,2,2)
        
        for n = 1:N
            plot(sqrtVAR(:,n),SimVarBLS(:,n),'s','Color',PlotOpts.colors(n,:))
            hold on
        end
        xlabel('Bias (ms)')
        ylabel('Expected bias (ms)')
        
        figure('Name','Actual vs averaging prediction')
        subplot(1,2,1)
        for n = 1:N
            plot(BIAS(:,n),SimBiasAve(:,n),'.','Color',PlotOpts.colors(n,:))
            hold on
        end
        xlabel('Bias (ms)')
        ylabel('Expected bias (ms)')
        
        subplot(1,2,2)
        
        for n = 1:N
            plot(sqrtVAR(:,n),SimVarAve(:,n),'s','Color',PlotOpts.colors(n,:))
            hold on
        end
        xlabel('Bias (ms)')
        ylabel('Expected bias (ms)')
        
            
        
    case {'No','no','n','N','NO'}
        
    otherwise
        error(['Plot option ' Plot ' not recognized!'])
end

