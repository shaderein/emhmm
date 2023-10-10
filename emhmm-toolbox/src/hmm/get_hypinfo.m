function hypinfo = get_hypinfo(learn_hyps, vbopt)
% get_hypinfo - get information about learnable hyperparameters
% [internal function]
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-08-04
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% 2018-02-12: added 'W0log' to use log(W0), and 'W0isqrt' to use sqrt(1/W0)
%             'W0' defaults to 'W0isqrt'

dim = length(vbopt.mu0);

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
      
    case 'beta0'
      hypinfo(k).hypname   = 'log(beta0)';
      hypinfo(k).optname   = 'beta0';
      hypinfo(k).derivname = 'd_logbeta0';
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
      
    case 'mu0'
      hypinfo(k).hypname   = 'm0';
      hypinfo(k).optname   = 'mu0';
      hypinfo(k).derivname = 'd_m0';
      hypinfo(k).hyptrans    = @(x) x;
      hypinfo(k).hypinvtrans = @(x) x;
      hypinfo(k).hypdims   = length(vbopt.mu0);
      
    otherwise
      error('bad value of learn_hyp');
      
  end
end
