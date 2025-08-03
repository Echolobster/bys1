function [w_int, RMSE, RMSEtest] = optimalW_INT(wm,wm_mem,samples,varargin)
%% optimalW_INT
%
%
%%

%% Defaults
% Estimator
estimator_default.type = 'BLS';

% Integration options
method_opts_default.dx = 10;
method_opts_default.type = 'quad';

% fminsearch options
OPTIONS_default = optimset('Display','iter');

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'wm')
addRequired(Parser,'wm_mem')
addRequired(Parser,'samples')
addParameter(Parser,'wp',0)
addParameter(Parser,'wVec',NaN)
addParameter(Parser,'estimator',estimator_default)
addParameter(Parser,'method_opts',method_opts_default)
addParameter(Parser,'trials',50000)
addParameter(Parser,'OPTIONS',OPTIONS_default)

parse(Parser,wm,wm_mem,samples,varargin{:})

wm = Parser.Results.wm;
wm_mem = Parser.Results.wm_mem;
samples = Parser.Results.samples;
wp = Parser.Results.wp;
wVec = Parser.Results.wVec;
estimator = Parser.Results.estimator;
method_opts = Parser.Results.method_opts;
trials = Parser.Results.trials;
OPTIONS = Parser.Results.OPTIONS;

%% Find optimal w_int
minimizant = @(p)(rmseW_INT(p,wm,wm_mem,wp,samples,trials,estimator,method_opts));

[w_int, RMSE] = fminsearch(minimizant,wm/sqrt(2),OPTIONS);

%% RMSE of w_ints
if isnan(wVec)
    wVec = linspace(0.01,max([wm wm_mem]),50);
end

RMSEtest = arrayfun(minimizant,wVec);


%% Functions

%% rmseW_INT
function rmse = rmseW_INT(w_int,wm,wm_mem,wp,samples,trials,estimator,...
    method_opts)

    % Set up estimator
    estimator.wm_drift = wm_mem;
    estimator.w_int = w_int;
    
    % Run simulation
    inds = ceil(length(samples)*rand(trials,1));
    for i = 1:length(inds)
        s(i,1) = samples(inds(i));
    end
    
    % Generate measuments of each ts
    noise = [wm_mem*s.*randn(size(s,1),1) wm*s.*randn(size(s,1),1)];
    m = repmat(s,1,2) + noise;
    
    % Generate estimates under model
    e(:,1) = ScalarBayesEstimators(m,wm,min(samples),max(samples),...
        'method',method_opts,'estimator',estimator);
    
    % Add production noise
    prodNoise = wp*e.*repmat(randn(size(e,1),1),1,size(e,2));
    p = e + prodNoise;
    
    % Find rmse
    rmse = sqrt( mean( (p-s).^2 ) );
    
    