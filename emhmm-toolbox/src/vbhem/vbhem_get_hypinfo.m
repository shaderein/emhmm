function hypinfo = vbhem_get_hypinfo(learn_hyps, vbopt)
% vbhem_get_hypinfo - get information about learnable hyperparameters
% [internal function]
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-07-29
% Antoni B. Chan, Janet H. Hsiao, Lan Hui
% City University of Hong Kong, University of Hong Kong

% 2021-07-29: v0.80 - initial version (adapted from vbhmm_get_hypinfo)

dim = length(vbopt.m0);

H = length(learn_hyps);
hypinfo = struct;
for k=1:H
  switch(learn_hyps{k})
    case 'alpha0'
      hypinfo(k).hypname   = 'log(alpha0)';
      hypinfo(k).optname   = 'alpha0';
      hypinfo(k).derivname = 'd_logalpha0';
      hypinfo(k).hyptrans    = @(x) exp(x);
      hypinfo(k).hypinvtrans = @(x) log(x);
      hypinfo(k).hypdims     = 1;
      
    case 'eta0'
      hypinfo(k).hypname   = 'log(eta0)';
      hypinfo(k).optname   = 'eta0';
      hypinfo(k).derivname = 'd_logeta0';
      hypinfo(k).hyptrans    = @(x) exp(x);
      hypinfo(k).hypinvtrans = @(x) log(x);
      hypinfo(k).hypdims     = 1;
      
    case 'epsilon0'
      hypinfo(k).hypname   = 'log(epsilon0)';
      hypinfo(k).optname   = 'epsilon0';
      hypinfo(k).derivname = 'd_logepsilon0';
      hypinfo(k).hyptrans    = @(x) exp(x);
      hypinfo(k).hypinvtrans = @(x) log(x);
      hypinfo(k).hypdims     = 1;
      
    case 'v0'
      hypinfo(k).hypname   = 'log(v0-D+1)';
      hypinfo(k).optname   = 'v0';
      hypinfo(k).derivname = 'd_logv0D1';
      hypinfo(k).hyptrans    = @(x) exp(x)+dim-1;
      hypinfo(k).hypinvtrans = @(x) log(x-dim+1);
      hypinfo(k).hypdims     = 1;
      
    case 'lambda0'
      hypinfo(k).hypname   = 'log(lambda0)';
      hypinfo(k).optname   = 'lambda0';
      hypinfo(k).derivname = 'd_loglambda0';
      hypinfo(k).hyptrans    = @(x) exp(x);
      hypinfo(k).hypinvtrans = @(x) log(x);
      hypinfo(k).hypdims     = 1;
    
    case {'W0', 'W0isqrt'}
      % sqrt transform
      hypinfo(k).hypname   = 'sqrt(W0inv)';
      hypinfo(k).optname   = 'W0';
      hypinfo(k).derivname = 'd_sqrtW0inv';
      hypinfo(k).hyptrans    = @(x) (x.^(-2));
      hypinfo(k).hypinvtrans = @(x) 1./sqrt(x);
      hypinfo(k).hypdims     = length(vbopt.W0);
    
    case 'W0log'
      % log transform
      hypinfo(k).hypname   = 'log(W0)';
      hypinfo(k).optname   = 'W0';
      hypinfo(k).derivname = 'd_logW0';
      hypinfo(k).hyptrans    = @(x) exp(x);
      hypinfo(k).hypinvtrans = @(x) log(x);
      hypinfo(k).hypdims     = length(vbopt.W0);
      
    case 'm0'
      hypinfo(k).hypname   = 'm0';
      hypinfo(k).optname   = 'm0';
      hypinfo(k).derivname = 'd_m0';
      hypinfo(k).hyptrans    = @(x) x;
      hypinfo(k).hypinvtrans = @(x) x;
      hypinfo(k).hypdims   = length(vbopt.m0);
      
    otherwise
      error('bad value of learn_hyp');
      
  end
end
