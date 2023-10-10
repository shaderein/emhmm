function [h3m_new] = form_outputH3M(h3m_r,Syn_STATS,Nj,hat_Z,L_elbo,covmode)
% (internal function)
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-07-29
% Antoni B. Chan, Janet H. Hsiao, Lan Hui
% City University of Hong Kong, University of Hong Kong

% 2021-07-29: v0.80 - initial version
%                   - remove L_elbo (takes too much space when replicated)

% using the hyperparameters to update paramters
[h3m_r] = convert_h3mrtoh3mb(h3m_r, covmode);

h3m_new = h3m_r;

h3m_new.Z = hat_Z;

h3m_new.hmm = {};
Kr = length(h3m_r.hmm);

for j = 1: Kr 
    
    % save variational parameters % for compute DIC
    % using the hmm hyperparamter name
    
    h3m_new.hmm{j}.varpar.alpha = h3m_r.hmm{j}.eta;
    h3m_new.hmm{j}.varpar.epsilon = h3m_r.hmm{j}.epsilon;
    
    Sr = length(h3m_r.hmm{j}.eta);
    
    for k =1:Sr
    h3m_new.hmm{j}.varpar.beta(k,1) = h3m_r.hmm{1, j}.emit{1, k}.lambda;
    h3m_new.hmm{j}.varpar.v(k,1) =h3m_r.hmm{1, j}.emit{1, k}.v;
    h3m_new.hmm{j}.varpar.m(:,k) = h3m_r.hmm{1, j}.emit{1, k}.m;
    h3m_new.hmm{j}.varpar.W(:,:,k) = h3m_r.hmm{1, j}.emit{1, k}.W;
    end
    

    h3m_new.hmm{j}.N1 = Syn_STATS{j}.Nj_rho1;
    h3m_new.hmm{j}.M = Syn_STATS{j}.Nj_rho2rho;
    h3m_new.hmm{j}.N = Syn_STATS{j}.Nj_rho;
       
    h3m_new.hmm{j}.prior = h3m_r.hmm{1, j}.prior;
    h3m_new.hmm{j}.trans = h3m_r.hmm{1, j}.A;
        
    h3m_new.hmm{j}.pdf = {};   
    for k = 1:Sr 
        h3m_new.hmm{j}.pdf{1,k}.mean = h3m_r.hmm{1, j}.emit{1, k}.centres;
        h3m_new.hmm{j}.pdf{1,k}.cov = h3m_r.hmm{1, j}.emit{1, k}.covars;
        h3m_new.hmm{j}.pdf{1,k}.covar_type = h3m_r.hmm{1, j}.emit{1, k}.covar_type;
    end
end
 
h3m_new.Nj = Nj;

%h3m_new.L_elbo = L_elbo;
 
% get cluster assignments
[foo, maxZ] = max(hat_Z, [], 2);
h3m_new.label = maxZ(:)';
 
% get cluster memberships
h3m_new.groups = {};
h3m_new.group_size = [];
for j=1:Kr
  h3m_new.groups{j} = find(h3m_new.label == j);
  h3m_new.group_size(j) = length(h3m_new.groups{j});
end

