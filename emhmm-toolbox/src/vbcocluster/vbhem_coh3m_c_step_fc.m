function [h3m_rs2dim] = vbhem_coh3m_c_step_fc(h3m_rs,h3m_bs,vbhemopt) 
% (internal function)
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-07-29
% Antoni B. Chan, Janet H. Hsiao, Lan Hui
% City University of Hong Kong, University of Hong Kong

% 2021-07-29: v0.80 - Hui - initial version


VERBOSE_MODE = vbhemopt.verbose;
 % start time
covmode = vbhemopt.emit.covar_type;

h3m_rs2dim.h3m_rs = h3m_rs;
h3m_rs2dim.K = h3m_rs(1).K;

% number of components in the base and reduced mixtures
Kb = vbhemopt.subjectsize;
Kr = h3m_rs2dim.K;
Nstimuli = vbhemopt.stimulisize; 

% number of states 
% Sbs = cell(Nstimuli,1);
% for i=1:Nstimuli
%     tmpSb =  cellfun(@(x) length(x.prior), h3m_bs(i).hmm);
%     Sbs{i} = tmpSb;
% end

if length(vbhemopt.S)>1
    Srs = vbhemopt.S;
else
    Srs = vbhemopt.S*ones(Nstimuli,1);
end

% length of the virtual sample sequences
T = vbhemopt.tau;
 
dim = h3m_bs(1).hmm{1}.emit{1}.nin;
% number of virtual samples 
virtualSamples = vbhemopt.Nv;

% correct
% ABC 2020-08-06: BUG FIX v0.xx: compute N_i for each stimuli
% it can be different because some trials may be missing for different stimuli
% and omega=0 for those trials
tilde_N_k_i = cell(Nstimuli,1);
for k=1:Nstimuli
  tmp = virtualSamples * h3m_bs(k).omega * Kb;
  tilde_N_k_i{k} = tmp(:);
end
 
maxIter = vbhemopt.max_iter; %maximum iterations allowed
minDiff = vbhemopt.minDiff; %termination criterion
 
if ~isfield(vbhemopt,'min_iter')
  vbhemopt.min_iter = 0;
end

% calculate derivative?
if vbhemopt.calc_LLderiv
  do_deriv = 1;
else
  do_deriv = 0;
end 
 
[vbhemopt, hyp_clipped] = vbhem_clip_hyps(vbhemopt);
vbhemopt.hyp_clipped = hyp_clipped;
    
% 2017-08-08 - ABC - handle MEX 
if ~isfield(vbhemopt, 'useMEX')
  vbhemopt.useMEX = 1;
end
% 
alpha0   = vbhemopt.alpha0;

canuseMEX = 0;
if (vbhemopt.useMEX)
  if (exist('vbhem_hmm_bwd_fwd_mex')==3)
    canuseMEX = 1;
    for i=1:length(h3m_bs)
        for n=1:length(h3m_bs(i).hmm)
          for j=1:length(h3m_bs(i).hmm{n}.emit)
            if h3m_bs(i).hmm{n}.emit{j}.ncentres ~= 1
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
    end
    
  else
    canuseMEX = 0;
    warning('vbhem_hmm_bwd_fwd_mex: MEX implementation not found');
  end
end


if (canuseMEX)
  % get maxN and maxN2
  maxNs = [];
  maxN2s = [];
  for i=1:vbhemopt.stimulisize
      maxN  = max(cellfun(@(x) length(x.prior), h3m_bs(i).hmm));
      maxN2 = max(cellfun(@(x) length(x.eta), h3m_rs(i).hmm));
      maxNs = [maxNs; maxN];
      maxN2s = [maxN2s; maxN2];
  end
  
else
  % 2017-08-07 - ABC - Caching/reshaping for faster E-step
  h3m_b_hmm_emit_cache = cell(emopt.stimulisize,Kb);
  for k=1:emopt.stimulisize
      for i=1:Kb
        g3m_b = h3m_bs(k).hmm{i}.emit;
        BKall = cellfun(@(x) x.ncentres, g3m_b);
        BKmax = max(BKall);

        % does not support different number of centres (for speed)
        if ~all(BKall==BKmax)
          error('different number of centres is not supported.');
        end

        gmmBcentres_tmp = cellfun(@(x) x.centres, g3m_b, ...
          'UniformOutput', false);
        h3m_b_hmm_emit_cache{k,i}.gmmBcentres = cat(3,gmmBcentres_tmp{:});  % BKmax x dim x N

        % extract all the covars
        gmmBcovars_tmp = cellfun(@(x) x.covars, g3m_b, ...
          'UniformOutput', false);
        switch (g3m_b{1}.covar_type)
          case 'diag'
            h3m_b_hmm_emit_cache{k,i}.gmmBcovars = cat(3,gmmBcovars_tmp{:});  % BKmax x dim x N
          case 'full'
            h3m_b_hmm_emit_cache{k,i}.gmmBcovars = cat(4,gmmBcovars_tmp{:});  % dim x dim x BKmax x N
            % TODO: add cache for logdet(gmmBcovars)
        end

        % extract all priors
        gmmBpriors_tmp = cellfun(@(x) x.priors', g3m_b, ...
          'UniformOutput', false);
        h3m_b_hmm_emit_cache{k,i}.gmmBpriors = cat(2,gmmBpriors_tmp{:});  % BKmax x N
      end
  end
end

% initialize arrays for E-step
L_elbo           = zeros(Nstimuli,Kb,Kr);
nu_1             = cell(Kb,Kr,Nstimuli);
update_emit_pr   = cell(Kb,Kr,Nstimuli);
update_emit_mu   = cell(Kb,Kr,Nstimuli);
update_emit_M    = cell(Kb,Kr,Nstimuli);
sum_xi           = cell(Kb,Kr,Nstimuli);

lastL = -realmax; %log-likelihood
 
iter = 0;
  
while 1
    %% E step
    %compute the responsibilities z_hat, phi_hat and statistics and constants
    
    if (canuseMEX)        
    %% MEX implementation   
    for s = 1 : Nstimuli
        
        Sr = Srs(s);
        for j = 1 : Kr
            
            hmm_r   = h3m_rs(s).hmm{j};
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
              %  m = hmm_r.emit{1, k}.m;
                
                t1 = psi(0, 0.5*repmat(v+1,dim,1) - 0.5*[1:dim]');
                
                switch(covmode)
                    case 'diag'
                        logLambdaTilde(k) = sum(t1) + const  + sum(log(diag(W)));

                    case 'full'
                        logLambdaTilde(k) = sum(t1) + const  + log(det(W));           
                end
                
                h3m_rs(s).hmm{j}.emit{k}.logLambdaTilde = logLambdaTilde(k);
                
                % -Elog(|Lambda|) + d/lambda
                % just here need
                h3m_rs(s).hmm{j}.emit{k}.logLambdaTildePlusDdivlamda =  -logLambdaTilde(k) + dim/lambda;
                
                psiEpsilonHat(k)  = psi(0,sum(epsilon(k,:)));
                logATilde(k,:)    = psi(0,epsilon(k,:)) - psiEpsilonHat(k);
            end
            
            psiEtaHat = psi(0,sum(eta));
            logPiTilde = psi(0,eta) - psiEtaHat;
            
            h3m_rs(s).hmm{j}.logPiTilde = logPiTilde;
            h3m_rs(s).hmm{j}.logATilde = logATilde;
        end
        
    end
        
        % 2017-11-25: ABC: BUG FIX: n->1, j->1.  just need to check covariance type
          switch(covmode)
              case 'diag'  
                  for k=1:Nstimuli
                      [L_elbo(k,:,:),nu_1(:,:,k), update_emit_pr(:,:,k), update_emit_mu(:,:,k), update_emit_M(:,:,k), sum_xi(:,:,k)] = ...
                          vbhem_hmm_bwd_fwd_mex(h3m_bs(k).hmm, h3m_rs(k).hmm,T,maxNs(k), maxN2s(k));
                  end
                  
              case 'full'
                  
                  for k=1:Nstimuli
                      
                      Sr = Srs(k);
                      logdetCovPlusDdivlamR = cell(1,Kr);
                      invCovR    = cell(1,Kr);
                      
                      for j=1:Kr
                          % does not support ncentres>1
                          logdetCovPlusDdivlamR{j} = zeros(1,length(h3m_rs(k).hmm{1, j}.eta));
                          invCovR{j}    = zeros(dim, dim, length(h3m_rs(k).hmm{1, j}.eta));
                          for n=1:Sr
                              logdetCovPlusDdivlamR{j}(n)  = h3m_rs(k).hmm{1, j}.emit{1, n}.logLambdaTildePlusDdivlamda;
                              invCovR{j}(:,:,n) = h3m_rs(k).hmm{1, j}.emit{1, n}.v.*h3m_rs(k).hmm{1, j}.emit{1, n}.W;
                          end
                      end
                      
                       [L_elbo(k,:,:), nu_1(:,:,k), update_emit_pr(:,:,k), update_emit_mu(:,:,k), update_emit_M(:,:,k), sum_xi(:,:,k)]  = ...
                          vbhem_hmm_bwd_fwd_mex(h3m_bs(k).hmm, h3m_rs(k).hmm,T,maxNs(k), maxN2s(k), logdetCovPlusDdivlamR, invCovR);
                  end

              otherwise  
                error('not supported')
          end
    
    else
        
        error('not support')
        %% compute phi(rho,beta)
        % loop reduced mixture components
%         
%         h3m_rs.hmms = h3m_rs.hmm;
%         
%         for j = 1 : Kr
%             
%             hmm_r   = h3m_rs.hmms{j};
%             Sr = Srs(j);
%             
%             eta     = hmm_r.eta;
%             epsilon = hmm_r.epsilon;
%             m       = hmm_r.m;
%             W       = hmm_r.W;
%             v       = hmm_r.v;
%             lambda  = hmm_r.lambda;
%             
%             const             = dim*log(2);
%             logLambdaTilde    = zeros(Sr,1);
%             logATilde         = zeros(Sr,Sr);
%             
%             for k = 1:Sr
%                 t1                = psi(0, 0.5*repmat(v(k)+1,dim,1) - 0.5*[1:dim]');
%                 logLambdaTilde(k) = sum(t1) + const  + log(det(W(:,:,k)));
%                 psiEpsilonHat(k)  = psi(0,sum(epsilon(k,:)));
%                 logATilde(k,:)    = psi(0,epsilon(k,:)) - psiEpsilonHat(k);
%             end
%             
%             psiEtaHat  = psi(0,sum(eta));
%             logPiTilde = psi(0,eta) - psiEtaHat;
%             
%             % loop base mixture components
%             for i = 1 : Kb
%                 
%                 hmm_b = h3m_bs.hmm{i};
%                 
%                 [L_elbo(i,j), bfopt{i,j}] = vbhem_hmm_bwd_fwd_fast(hmm_b,T,m,...
%                     W,v,lambda,logLambdaTilde,logATilde,logPiTilde, h3m_b_hmm_emit_cache{i});
%                 % end loop base
%             end
%             % end loop reduced
%             
%             % just save in struct of i=1
%             bfopt{1,j}.logLambdaTilde = logLambdaTilde;
%             bfopt{1,j}.logATilde = logATilde;
%             bfopt{1,j}.logPiTilde = logPiTilde;
%         end
        
    end % end canuseMEX
  
    %% compute hat_z
    %logOmegaTilde = zeros(Nstimuli,Kr);
    % the first h3m_rs(i).alpha (initialization), for each i, it is different
    % if initialization is vbh3m, it will be the same
    alpha = zeros(1,Kr);
    for i =1: Nstimuli
        alpha = alpha + h3m_rs(i).alpha;
    end
    
    alpha = alpha/Nstimuli;
    psiAlphaHat = psi(0,sum(alpha));
    logOmegaTilde = psi(0,alpha) - psiAlphaHat;
 
    % compute the z_ij  
    log_Z = zeros(Kb,Kr);
    sum_s_Ni = zeros(Kb,1);
    for k=1:Nstimuli
        log_Z = log_Z +  tilde_N_k_i{k}(:).* reshape(L_elbo(k,:,:),[Kb,Kr]);
        sum_s_Ni = sum_s_Ni + tilde_N_k_i{k}(:);
    end
    Sum_s_L_elbo = log_Z;
    
    log_Z = log_Z + (sum_s_Ni * ones(1,Kr)) .* (ones(Kb,1) * logOmegaTilde);
  
    hat_Z = exp(log_Z - logtrick(log_Z')' * ones(1,Kr)); 
    hat_Z = hat_Z + 1e-50;
    
    % scale the z_ij by the number of virtual samples
    Z = zeros(Kb,Kr,Nstimuli);
    for s = 1: Nstimuli
        Z(:,:,s) = hat_Z .* (tilde_N_k_i{s}(:)* ones(1,Kr));%(sum_s_Ni * ones(1,Kr));
    end
    
    NL = sum(sum(Z,3),1);       %  sum(hat_Z,1)
    NL = NL + 1e-50;
    
    %% calculate the variational lower bound %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
    h3m_stats               ={};
    h3m_stats.Srs           = Srs;
    h3m_stats.logOmegaTilde = logOmegaTilde;
    h3m_stats.hat_Z         = hat_Z;
    h3m_stats.NL            = NL;
    h3m_stats.do_deriv      = do_deriv;
    h3m_stats.Sum_s_L_elbo  = Sum_s_L_elbo;

 
    if (canuseMEX)  
       [L] = vbhemh3m_colb(h3m_stats, h3m_rs, vbhemopt);       
    else
      error('not support') %[L,~,Lsave] = vbhemh3m_lb(h3m_stats, h3m_rs, vbhemopt, bfopt);
    end
       
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
      if canuseMEX      
        [Lnew, dL] = vbhemh3m_colb(h3m_stats, h3m_rs, vbhemopt) ;
      else
        [Lnew, dL] = vbhemh3m_colb(h3m_stats, h3m_rs, vbhemopt, bfopt) ;
      end
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

   %% M-step
   % update hyperparameters
   Syn_STATS = cell(Kr,Nstimuli);
   
   for s = 1:Nstimuli
       
    [h3m_r] = convert_h3mrtoh3mb(h3m_rs(s), covmode);  
    alpha = alpha0 + NL;
    h3m_r.alpha = alpha;
    
    for j = 1:Kr
        
        estats.Z_Ni    = Z(:,j,s);
        estats.nu_1    = nu_1(:,j,s);%nu_1(:,j);
        estats.sum_xi  = sum_xi(:,j,s);%sum_xi(:,j);
        estats.emit_pr = update_emit_pr(:,j,s);
        estats.emit_mu = update_emit_mu(:,j,s);
        estats.emit_M  = update_emit_M(:,j,s);
        
        % setup constants
        estats.Sr  = Srs(s);
        estats.dim = dim;
        estats.Kb  = Kb;
        
        if ~do_break            
            [h3m_r.hmm{j}]= vbhem_mstep_component(estats, vbhemopt);
        else
            [h3m_r.hmm{j},Syn_STATS{j,s}]= vbhem_mstep_component(estats, vbhemopt);
        end        
    end

    h3m_rs(s).hmm =h3m_r.hmm;
    h3m_rs(s).alpha = h3m_r.alpha;
    
   end
     
   %% check convergency
     
    iter = iter +1;
    lastL = L;  
 
   % break out of E & M loop
    if (do_break)
      break;
    end
    
end

%% generate the output h3m
h3m_new     = cell(1,Nstimuli);

for s =1:Nstimuli   
    [h3m_new{s}] = form_outputH3M(h3m_rs(s),Syn_STATS(:,s),NL,hat_Z,L_elbo,covmode);
end
% 
clear h3m_rs

h3m_rs2dim.h3m_rs = h3m_new;
h3m_rs2dim.label = h3m_new{1, 1}.label;
h3m_rs2dim.groups = h3m_new{1, 1}.groups;
h3m_rs2dim.group_size = h3m_new{1, 1}.group_size;

h3m_rs2dim.L_elbo = L_elbo;

 % save derivatives
if (do_deriv)
  h3m_rs2dim.dLL = dL;
end 
 
h3m_rs2dim.iter = iter;
h3m_rs2dim.LL = lastL;
h3m_rs2dim.NL = NL;
