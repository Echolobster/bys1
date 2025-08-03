function BIC = ModelBIC(varargin)
%% ModelBIC
%
%   BIC = ModelBIC()
%
%   Evaluates the Bayesian information criterion for a data set given a set
%   of model parameters fit to data.
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'param',NaN)

parse(Parser,varargin{:})

param = Parser.Results.param;

%% Evaluate log likelihood of the data, given model parameters
switch model
    case 'BLSbiasedLapse'
        llikelihood = BLSbiasedLapse_Validator(s,m,wm,wp,b,lapse
        
    case 'SubOptMemBiasbiasedLapse'
        
    case 'MAPbiasedLapse'
        
    otherwise
        error(['Bayesian information criterion not yet supported for model type ' model '!'])
        
end