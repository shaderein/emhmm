function [h3m_out] = vbhem_coh3m_c1(h3m_bs, vbhemopt, K) 
% (internal function)
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-07-29
% Antoni B. Chan, Janet H. Hsiao, Lan Hui
% City University of Hong Kong, University of Hong Kong

% 2021-07-29: v0.80 - Hui - initial version

vbhemopt.K = K;
%vbhemopt.S = S;


if ~strcmp(vbhemopt.initmode, 'auto')
    
    %% run VBHEM clustering

    [h3m_out] = vbhem_coh3m_c(h3m_bs, vbhemopt);
          
else
    %% auto initialization - try these methods

    if isfield(vbhemopt,'initmodes')
        initmodes = vbhemopt.initmodes;
        initopts = vbhemopt.initopts;
    else
        initmodes = { 'baseem', 'gmmNew', 'gmmNew2','gmmNew2_r','wtkmeans'};
        initopts = { 'u', 'r0', 'u0','u0','r0'};
    end
    
    %          initmodes = {'random', 'baseem', 'gmmNew', 'gmmNew2','gmmNew2_r'};
    %          initopts = {'u', 'u', 'r0', 'u0','u0'};
    % v0.80 - ABC - use normal for loop. parfor happens in vbhem_coh3m_c
    for ii=1:length(initmodes)
        % tic;
        vbhemopt3 = vbhemopt;
        vbhemopt3.initmode = initmodes{ii};
        vbhemopt3.initopt.mode = initopts{ii};
        fprintf('auto initialization: trying %s: ', vbhemopt3.initmode);
        h3m_out_trials{ii} = vbhem_coh3m_c(h3m_bs, vbhemopt3);
        h3m_out_trials{ii}.Initmodes = initmodes{ii};
        trials_LL(ii) = h3m_out_trials{ii}.LL;
    %    h3m_out_trials{ii}.time = toc;
    end
    % get best method
    [maxLL, ind] = max(trials_LL);
    h3m_out = h3m_out_trials{ind};
    h3m_out.h3m_out_trials = h3m_out_trials;
    fprintf('best init was %s; LL=%g\n', initmodes{ind}, trials_LL(ind));
end