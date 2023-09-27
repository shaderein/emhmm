function d = computeCohen_d(x1, x2, varargin)
% 
% call: d = computeCohen_d(x1, x2, varargin)
% 
% EFFECT SIZE of the difference between the two 
% means of two samples, x1 and x2 (that are vectors), 
% computed as "Cohen's d". 
% 
% If x1 and x2 can be either two independent or paired 
% samples, and should be treated accordingly:
%  
%   d = computeCohen_d(x1, x2, 'independent');  [default]
%   d = computeCohen_d(x1, x2, 'paired');
% 
% Note: according to Cohen and Sawilowsky:
%
%      d = 0.01  --> very small effect size
%      d = 0.20  --> small effect size
%      d = 0.50  --> medium effect size
%      d = 0.80  --> large effect size
%      d = 1.20  --> very large effect size
%      d = 2.00  --> huge effect size
%
%
% Ruggero G. Bettinardi (RGB)
% Cellular & System Neurobiology, CRG
% -------------------------------------------------------------------------------------------
%
% Code History:
%
% 25 Jan 2017, RGB: Function is created

% - v0.78: initial version
  
if nargin < 3, testType = 'independent'; 
else           testType = varargin{1}; 
end
% basic quantities:
n1       = numel(x1);
n2       = numel(x2);
mean_x1  = nanmean(x1);
mean_x2  = nanmean(x2);
var_x1   = nanvar(x1);
var_x2   = nanvar(x2);
meanDiff = (mean_x1 - mean_x2);
% select type of test:
isIndependent = strcmp(testType, 'independent');
isPaired      = strcmp(testType, 'paired');
% compute 'd' accordingly:
if isIndependent
    
    sv1      = ((n1-1)*var_x1);
    sv2      = ((n2-1)*var_x2);
    numer    =  sv1 + sv2;
    denom    = (n1 + n2 - 2);
    pooledSD =  sqrt(numer / denom); % pooled Standard Deviation
    s        = pooledSD;             % re-name
    d        =  meanDiff / s;        % Cohen's d (for independent samples)
    
elseif isPaired
    
    haveNotSameLength = ~isequal( numel(x1), numel(x2) );
    if haveNotSameLength, error('In a paired test, x1 and x2 have to be of same length!'), end
    
    deltas   = x1 - x2;         % differences
    sdDeltas = nanstd(deltas);  % standard deviation of the diffferences
    s        = sdDeltas;        % re-name
    d        =  meanDiff / s;   % Cohen's d (paired version)
    
end