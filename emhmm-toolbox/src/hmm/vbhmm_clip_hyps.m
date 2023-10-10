function [vbopt, clipped] = vbhmm_clip_hyps(vbopt)
% vbhmm_clip_hyps - clip hyperparameters at maximum/minimum levels for stable computations
%                   when doing hyperparameter estimation.
%
% this is an internal function called by vbhmm_em
%
% [vbopt, clipped] = vbhmm_clip_hyps(vbopt)
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-08-21
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% 2018-02-11: v0.72 added W0 (and diag W0)
%           : v0.72 use vbopt structure for min/max values

dim = length(vbopt.mu0);

clipped = struct;

%% REMINDER: change code in vbhmm_em_lb for clipping the derivatives!

% get from vbopt
MAX = vbopt.hyps_max;
MIN = vbopt.hyps_min;

% maximum values of hyperparameters, for stability
%MAX.alpha0   = 1.0686e+13;  % exp(30)
%MAX.epsilon0 = 1.0686e+13;  % exp(30)
%MAX.v0       = 1e4;
%MAX.beta0    = 1.0686e+13;  % exp(30)
%MAX.W0       = 1e13;  

% minimum values
%MIN.alpha0   = 1.0686e-13;  % exp(-30)
%MIN.epsilon0 = 1.0686e-13;  % exp(-30)
%MIN.v0       = 2.0612e-09+dim-1;  % exp(-20)+dim-1 - should be okay to ~10000 dimensions
%MIN.beta0    = 1.0686e-13;  % exp(-30)
%MIN.W0       = 1.0686e-13;

allhyps  = fieldnames(MAX);
allhyps2 = fieldnames(MIN);

% initialize clipped to 0
for h=1:length(allhyps)
  myhyp = allhyps{h};
  clipped.(myhyp) = zeros(1,length(vbopt.(myhyp)));
end
for h=1:length(allhyps2)
  myhyp = allhyps2{h};
  clipped.(myhyp) = zeros(1,length(vbopt.(myhyp)));
end

% check for max clip
for h=1:length(allhyps)
  myhyp = allhyps{h};  
  for i=1:length(vbopt.(myhyp))
    if (vbopt.(myhyp)(i) >= MAX.(myhyp))
      vbopt.(myhyp)(i) = MAX.(myhyp);
      clipped.(myhyp)(i) = +1;
      if (vbopt.verbose>1)
        fprintf('[MAX clipped %s(%d) to %g]', myhyp, i, MAX.(myhyp));
      end
    end
  end
end

% check for min clip
for h=1:length(allhyps2)
  myhyp = allhyps2{h};  
  for i=1:length(vbopt.(myhyp))
    if (vbopt.(myhyp)(i) <= MIN.(myhyp))
      vbopt.(myhyp)(i) = MIN.(myhyp);
      clipped.(myhyp)(i) = -1;
      if (vbopt.verbose>1)
        fprintf('[MIN clipped %s(%d) to %g]', myhyp, i, MIN.(myhyp));
      end
    end
  end  
end

