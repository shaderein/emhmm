function [pyy] = vbhmm_prob_jointobs2(hmm, t)
% vbhmm_prob_jointobs2 - compute joint distribution of two observations in an HMM
%
%   [pxx] = vbhmm_prob_jointobs2(hmm, t)
%
%   Compute the joint distribution p(x_{t-1}, x_t)
%
%   INPUT:  hmm - the hmm
%             t - the time step (fixation number). 
%                 inf = use steady-state.
%
%   OUTPUT: pxx - p(x_{t-1}, x_t)
%                 (a 4-dim GMM, where the variable is [x_{t-1}, x{t}]
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2018-10-24
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% VERSIONS
% 2021-07-21: v0.78 ABC - initial version

if t<2
  error('t must be >= 2');
end

S = length(hmm.pdf);
D = length(hmm.pdf{1}.mean);

if isinf(t)
  % get the steady state probabilites
  pz = vbhmm_prob_steadystate(hmm);
else
  % get state marginal at t-1
  tmp = vbhmm_prob_marginalstate(hmm, t-1);
  pz = tmp(:,end);
end
  
% get the joint state probability
% p(z_{t-1}, z_t)
pzz = hmm.trans .* repmat(pz(:), [1 S]);

% build p(y_{t-1}, y_t)
pyy.prior = zeros(1,S*S);
pyy.pdf   = cell(1,S*S);

ind = 1;
for i=1:S    % t-1
  for j=1:S  % t
    tmpg.mean = [hmm.pdf{i}.mean, hmm.pdf{j}.mean];
    tmpg.cov  = [hmm.pdf{i}.cov, zeros(D); zeros(D), hmm.pdf{j}.cov];    
    pyy.prior(ind) = pzz(i,j);
    pyy.pdf{ind}   = tmpg;    
    ind = ind+1;
  end
end



