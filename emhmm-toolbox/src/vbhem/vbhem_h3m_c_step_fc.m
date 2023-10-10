function [h3m_r_new] = vbhem_h3m_c_step_fc(h3m_r,h3m_b,vbhemopt) 
% (internal function)
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-07-29
% Antoni B. Chan, Janet H. Hsiao, Lan Hui
% City University of Hong Kong, University of Hong Kong

% 2021-07-29: v0.80 - initial version

% if 'diag', W also is 'diag' [1,dim], m: [1,dim]
% if 'full', W also is 'full' [dim,dim]

VERBOSE_MODE = vbhemopt.verbose;

covmode = vbhemopt.emit.covar_type;
h3m_r.LogLs = [];

% number of components in the base and reduced mixtures
Kb = h3m_b.K;
Kr = h3m_r.K;
 
% number of states 
Sr = h3m_r.S;  
 
% length of the virtual sample sequences
T = vbhemopt.tau;
 
% dimension of the emission variable
dim = h3m_b.hmm{1}.emit{1}.nin;
 
% number of virtual samples 
virtualSamples = vbhemopt.Nv  * Kb;

% correct
tilde_N_k = virtualSamples * h3m_b.omega; 
tilde_N_k = tilde_N_k(:);%%%change n_i to lie vector;
 
maxIter = vbhemopt.max_iter; %maximum iterations allowed
minDiff = vbhemopt.minDiff; %termination criterion
 
lastL = -realmax; %log-likelihood
 
iter = 0;
 
if ~isfield(vbhemopt,'min_iter')
  vbhemopt.min_iter = 0;
end

if vbhemopt.calc_LLderiv
  do_deriv = 1;
else
  do_deriv = 0;
end 
    
% 2017-08-08 - ABC - handle MEX 
if ~isfield(vbhemopt, 'useMEX')
  vbhemopt.useMEX = 1;
end
 
alpha0   = vbhemopt.alpha0;      

canuseMEX = 0;
if (vbhemopt.useMEX)
  if (exist('vbhem_hmm_bwd_fwd_mex')==3)
    canuseMEX = 1;
    for n=1:length(h3m_b.hmm)
      for j=1:length(h3m_b.hmm{n}.emit)
        if h3m_b.hmm{n}.emit{j}.ncentres ~= 1
          canuseMEX = 0;
          warning('vbhem_hmm_bwd_fwd_mex: ncentres should be 1');
        end

        if canuseMEX == 0
          break;
        end
      end
      if canuseMEX == 0
        break
      end
    end
    
  else
    canuseMEX = 0;
    warning('vbhem_hmm_bwd_fwd_mex: MEX implementation not found');
  end
end

if (canuseMEX)
  % get maxN and maxN2
  maxN  = max(cellfun(@(x) length(x.prior), h3m_b.hmm));
  maxN2 = max(cellfun(@(x) length(x.eta), h3m_r.hmm));  %eta
  
else
  
  % make cache for faster E-step
  % v0.73 - put into its own function
  h3m_b_hmm_emit_cache = make_h3m_b_cache(h3m_b);
  
  % initialize arrays for E-step 
  % start looping variational E step and M step
%   L_elbo  = zeros(Kb,Kr);
%   nu_1             = cell(Kb,Kr);
%   update_emit_pr   = cell(Kb,Kr);
%   update_emit_mu   = cell(Kb,Kr);
%   update_emit_M    = cell(Kb,Kr);
%   sum_xi           = cell(Kb,Kr);
  
  xL_elbo  = zeros(Kb,Kr);
  xnu_1             = cell(Kb,Kr);
  xupdate_emit_pr   = cell(Kb,Kr);
  xupdate_emit_mu   = cell(Kb,Kr);
  xupdate_emit_M    = cell(Kb,Kr);
  xsum_xi           = cell(Kb,Kr);
  
  
end

while 1
    %% E step
    %compute the responsibilities z_hat, phi_hat and statistics and constants
    for j = 1 : Kr
        
        hmm_r   = h3m_r.hmm{j};
        epsilon = hmm_r.epsilon;
        eta     = hmm_r.eta; 
        
        const             = dim*log(2);      
        logLambdaTilde    = zeros(Sr,1);
        logATilde         = zeros(Sr,Sr);
        psiEpsilonHat     = zeros(1,Sr);
        
        for k = 1:Sr
            
            v = hmm_r.emit{k}.v;
            lambda = hmm_r.emit{1, k}.lambda;
            W = hmm_r.emit{1, k}.W;
            m = hmm_r.emit{1, k}.m;
            
            t1 = psi(0, 0.5*repmat(v+1,dim,1) - 0.5*[1:dim]');            
            
            switch(covmode)
                case 'diag'
                   logLambdaTilde(k) = sum(t1) + const  + sum(log(diag(W)));

                   % h3m_r.hmm{j}.emit{k}.covars = diag(Cov_new)';
                    
                case 'full'
                   % h3m_r.hmm{j}.emit{k}.covars = Cov_new;
                    logLambdaTilde(k) = sum(t1) + const  + log(det(W));

            end
            
            h3m_r.hmm{j}.emit{k}.logLambdaTilde = logLambdaTilde(k);
            
            % -Elog(|Lambda|) + d/lambda
            % just here need
            h3m_r.hmm{j}.emit{k}.logLambdaTildePlusDdivlamda =  -logLambdaTilde(k) + dim/lambda;

            psiEpsilonHat(k)  = psi(0,sum(epsilon(k,:)));
            logATilde(k,:)    = psi(0,epsilon(k,:)) - psiEpsilonHat(k);
        end
        
        psiEtaHat = psi(0,sum(eta));
        logPiTilde = psi(0,eta) - psiEtaHat;
        
        h3m_r.hmm{j}.logPiTilde = logPiTilde;
        h3m_r.hmm{j}.logATilde = logATilde;     
    end
    
    
    if (canuseMEX)
        %% MEX implementation
        
        % 2017-11-25: ABC: BUG FIX: n->1, j->1.  just need to check covariance type
        switch(h3m_b.hmm{1}.emit{1}.covar_type)
            case 'diag'
               
                [L_elbo, nu_1, update_emit_pr, update_emit_mu, update_emit_M, sum_xi] = ...
                    vbhem_hmm_bwd_fwd_mex(h3m_b.hmm,h3m_r.hmm,T,maxN, maxN2);
                
            case 'full'
                % cache logdets and inverse covariances
                logdetCovPlusDdivlamR = cell(1,Kr);
                invCovR    = cell(1,Kr);
                
                for j=1:Kr
                    % does not support ncentres>1
                    logdetCovPlusDdivlamR{j} = zeros(1,length(h3m_r.hmm{1, j}.eta));
                    invCovR{j}    = zeros(dim, dim, length(h3m_r.hmm{1, j}.eta));
                    for n=1:length(h3m_r.hmm{1, j}.eta)
                        logdetCovPlusDdivlamR{j}(n)  = h3m_r.hmm{1, j}.emit{1, n}.logLambdaTildePlusDdivlamda;
                        invCovR{j}(:,:,n) = h3m_r.hmm{1, j}.emit{1, n}.v.*h3m_r.hmm{1, j}.emit{1, n}.W;
                    end
                end
                
                [L_elbo, nu_1, update_emit_pr, update_emit_mu, update_emit_M, sum_xi] = ...
                    vbhem_hmm_bwd_fwd_mex(h3m_b.hmm,h3m_r.hmm,T,maxN, maxN2, logdetCovPlusDdivlamR, invCovR);
                
            otherwise
                error('not supported')
        end
        
    else
        % cache the variational parameters of h3m_r
        %% compute phi(rho,beta)
        % loop reduced mixture components    
        for j = 1 : Kr            
            hmm_r   = h3m_r.hmm{j};
            m       = zeros(dim,Sr);
            W       = zeros(dim,dim,Sr);
            v       = zeros(Sr,1);
            lambda  = zeros(Sr,1);
            logLambdaTilde    = zeros(Sr,1);          
            
            for k = 1:Sr
                m(:,k) = hmm_r.emit{1, k}.m';
                v(k) = hmm_r.emit{1, k}.v;
                switch(covmode)
                    case 'diag'
                        W(:,:,k) = diag(hmm_r.emit{1, k}.W);                                                
                    case 'full'
                        W(:,:,k) = hmm_r.emit{1, k}.W;       
                end
                
                lambda(k) = hmm_r.emit{1, k}.lambda;
                logLambdaTilde(k) = h3m_r.hmm{j}.emit{k}.logLambdaTilde;
            end
            
            logATilde   = hmm_r.logATilde;
            logPiTilde = hmm_r.logPiTilde;
            
            % loop base mixture components
            for i = 1 : Kb
                
                hmm_b = h3m_b.hmm{i};
                %                 [L_elbo(i,j) ...
                %                     nu_1{i,j} ...
                %                     update_emit_pr{i,j} ...
                %                     update_emit_mu{i,j} ...
                %                     update_emit_M{i,j}  ...
                %                     sum_xi{i,j} ...
                %                  ] = vbhem_hmm_bwd_fwd_fast(hmm_b,T,m,W,v,lambda,...
                %                  logLambdaTilde,logATilde,logPiTilde, h3m_b_hmm_emit_cache{i});
                [xL_elbo(i,j) ...
                    xnu_1{i,j} ...
                    xupdate_emit_pr{i,j} ...
                    xupdate_emit_mu{i,j} ...
                    xupdate_emit_M{i,j}  ...
                    xsum_xi{i,j} ...
                    ] = vbhem_hmm_bwd_fwd_fast(hmm_b,T,m,W,v,lambda,...
                    logLambdaTilde,logATilde,logPiTilde, h3m_b_hmm_emit_cache{i});
                
                %{
                err = [];
                err(end+1) = abs(xL_elbo(i,j)-L_elbo(i,j));
                err(end+1) = sum(abs(xnu_1{i,j}(:) - nu_1{i,j}(:)));
                err(end+1) = sum(abs(xupdate_emit_pr{i,j}(:) - update_emit_pr{i,j}(:)));
                err(end+1) = sum(abs(xupdate_emit_mu{i,j}(:) - update_emit_mu{i,j}(:)));
                err(end+1) = sum(abs(xupdate_emit_M{i,j}(:)  - update_emit_M{i,j}(:)));
                err(end+1) = sum(abs(xsum_xi{i,j}(:) - sum_xi{i,j}(:)));
                if any(err>1e-6)
                    warning('mismatch')
                    keyboard;
                end
                %}
                
                
                % end loop base
            end
        end  % end loop reduced
        
    end % end canuseMEX
    
    %% compute hat_z
    alpha = h3m_r.alpha;
    psiAlphaHat = psi(0,sum(alpha));
    logOmegaTilde = psi(0,alpha) - psiAlphaHat;
 
    
    log_Z = (tilde_N_k(:) * ones(1,Kr)).*(ones(Kb,1)* logOmegaTilde + L_elbo);
    hat_Z = exp(log_Z - logtrick(log_Z')' * ones(1,Kr)); 
    hat_Z = hat_Z + 1e-50;
    
    % scale the z_ij by the number of virtual samples
    Z_Ni = hat_Z .* (tilde_N_k(:) * ones(1,Kr));

    Nj = sum(Z_Ni,1);       
    Nj = Nj + 1e-50;    

     %% calculate the variational lower bound %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
    h3m_stats={};
    h3m_stats.logOmegaTilde = logOmegaTilde;
    h3m_stats.hat_Z         = hat_Z;
    h3m_stats.Z             = Z_Ni;
    h3m_stats.Nj            = Nj;
    h3m_stats.L_elbo        = L_elbo;
    h3m_stats.do_deriv      = do_deriv;
 
   
    [L] = vbhemh3m_lb(h3m_stats, h3m_r, vbhemopt);

    
    h3m_r.LL = L;  
    
    do_break = 0;
    unstable = 0;
    
    if iter == 1
      % show the first iteration too
      if (VERBOSE_MODE >= 3)
        fprintf('iter %d: L=%g; dL=x\n', 1, L);
      end
    end
    
    if iter > 1
      likIncr = abs((L-lastL)/lastL);
      if (VERBOSE_MODE >= 3)        
        fprintf('iter %d: L=%g; dL=%g', iter, L, likIncr);
        if (L-lastL < 0)
          fprintf('[LL decreased: %g]!!!', L-lastL);
          %keyboard
        end
        fprintf('\n');
      else 
        if (L-lastL < 0)
          if (VERBOSE_MODE >= 2)
            fprintf('[LL decreased]');
          end
        end
      end
      if (likIncr <= minDiff)
        do_break = 1;
      end
    end
    if (iter == maxIter)
      warning('max iterations reached');
      do_break = 1;
    end
    
    % if L is NaN, then usually it means an impossible model L=-inf.
  
    if any(isnan(L))
      if (VERBOSE_MODE >=1)
        warning('NaN found, unstable model');
      end
      do_break = 1;
      unstable = 1;
      
      % 2018-02-01: stability fix when optimizing hyps (when it tries weird parameter combos)
      L = -inf;
      
      % save for debugging
      if (VERBOSE_MODE >= 3)
        foo=tempname('.');
        warning('degenerate models detected: saving to %s', foo);
        save([foo '.mat'])
      end
    end
 
    % calculate derivatives if breaking (do this before the last M-step)
    if (do_break) && (do_deriv)
      h3m_stats.do_deriv = 1;
      [Lnew, dL] = vbhemh3m_lb(h3m_stats, h3m_r, vbhemopt) ;

      % invalidate the gradient
      if (unstable)
        tmp = fieldnames(dL);
        for i=1:length(tmp)
          dL.(tmp{i}) = nan*dL.(tmp{i});
        end
      end
    end
    
 
    % break early, since unstable model
    if (do_break) && unstable
      break;
    end
    
    %% compute the second phi
        
    % convert the h3m_r to the same form of h3m_b
    % at he same time, update the prior & A, center & coveras
    % SINCE we want to use hem_bwd_fwd
    [h3m_r] = convert_h3mrtoh3mb(h3m_r, covmode);
    
 %   [xL_elbo,nu_1,update_emit_pr,update_emit_mu,update_emit_M,sum_xi] = hem_bwd_fwd(h3m_b,h3m_r,T,vbhemopt);
%     
%     if iter >1
%         % xLelbo is bigger than Lelbo,
%         % L_elbo is a lower bound on L_elbo
%         check = sum(sum(L_elbo > xL_elbo))>1;
%         if check >1
%         keyboard
%         end
%     end
%     
    %% M-step
    % update hyperparameters
    alpha = alpha0 + Nj;
    h3m_r.alpha = alpha;
    Syn_STATS = cell(Kr,1);
    
    for j = 1:Kr
        
        estats.Z_Ni      = Z_Ni(:,j);
        estats.nu_1    = nu_1(:,j);
        estats.sum_xi  = sum_xi(:,j);
        estats.emit_pr = update_emit_pr(:,j);
        estats.emit_mu = update_emit_mu(:,j);
        estats.emit_M  = update_emit_M(:,j);
        
        % setup constants
        estats.Sr  =length(h3m_r.hmm{1, j}.eta);
        estats.dim = dim;
        estats.Kb  = Kb;
        
        if ~do_break            
            [h3m_r.hmm{j}]= vbhem_mstep_component(estats, vbhemopt);
        else
            [h3m_r.hmm{j},Syn_STATS{j}]= vbhem_mstep_component(estats, vbhemopt);
        end        
    end
    
 
   %% check convergency
     
    iter = iter +1;
    h3m_r.LogLs(end+1) = L;
    lastL = L;  
 
   % break out of E & M loop
    if (do_break)
      break;
    end
    
end
%% generate the output h3m
% update parameters
 [h3m_new] = form_outputH3M(h3m_r,Syn_STATS,Nj,hat_Z,L_elbo,covmode);
 
% save derivatives
if (do_deriv)
  h3m_new.dLL = dL;
end

% if (SAVE_PHI)
%     h3m_new.Theta_all = xTheta_all;
% end
% output
h3m_r_new = h3m_new;
    
    
    
    
    
    
 
 
 

