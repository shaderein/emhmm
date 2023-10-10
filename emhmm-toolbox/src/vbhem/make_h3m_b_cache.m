function h3m_b_hmm_emit_cache = make_h3m_b_cache(h3m_b)
% (internal function)
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-07-29
% Antoni B. Chan, Janet H. Hsiao, Lan Hui
% City University of Hong Kong, University of Hong Kong

% 2021-07-29: v0.80 - initial version

Kb = h3m_b.K;

% 2017-08-07 - ABC - Caching/reshaping for faster E-step
h3m_b_hmm_emit_cache = cell(1,Kb);
for i=1:Kb
  g3m_b = h3m_b.hmm{i}.emit;
  BKall = cellfun(@(x) x.ncentres, g3m_b);
  BKmax = max(BKall);
  
  % does not support different number of centres (for speed)
  if ~all(BKall==BKmax)
    error('different number of centres is not supported.');
  end
  
  gmmBcentres_tmp = cellfun(@(x) x.centres, g3m_b, ...
    'UniformOutput', false);
  h3m_b_hmm_emit_cache{i}.gmmBcentres = cat(3,gmmBcentres_tmp{:});  % BKmax x dim x N
  
  % extract all the covars
  gmmBcovars_tmp = cellfun(@(x) x.covars, g3m_b, ...
    'UniformOutput', false);
  switch (g3m_b{1}.covar_type)
    case 'diag'
      h3m_b_hmm_emit_cache{i}.gmmBcovars = cat(3,gmmBcovars_tmp{:});  % BKmax x dim x N
    case 'full'
      h3m_b_hmm_emit_cache{i}.gmmBcovars = cat(4,gmmBcovars_tmp{:});  % dim x dim x BKmax x N
      % TODO: add cache for logdet(gmmBcovars)
  end
  
  % extract all priors
  gmmBpriors_tmp = cellfun(@(x) x.priors', g3m_b, ...
    'UniformOutput', false);
  h3m_b_hmm_emit_cache{i}.gmmBpriors = cat(2,gmmBpriors_tmp{:});  % BKmax x N
end

end
