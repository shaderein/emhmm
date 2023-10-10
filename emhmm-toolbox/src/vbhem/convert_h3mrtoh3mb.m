function [h3m_new] = convert_h3mrtoh3mb(h3m_r, covmode)
% (internal function)
% convert the h3m_r to the same form of h3m_b
% covmode = {'diag','full'}
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-07-29
% Antoni B. Chan, Janet H. Hsiao, Lan Hui
% City University of Hong Kong, University of Hong Kong

% 2021-07-29: v0.80 - initial version

 h3m_new = h3m_r;
 
 K = length(h3m_r.hmm);
 dim = length(h3m_r.hmm{1, 1}.emit{1, 1}.m);
 
 for j = 1:K
    
    tempHmm = h3m_r.hmm{j};
    % prior
    eta_i = h3m_r.hmm{1, j}.eta;
    prior_tmp = sum( eta_i);
    prior = eta_i ./ prior_tmp;
    tempHmm.prior = prior';
    
    % A
    trans_t =  h3m_r.hmm{j}.epsilon;  % [row]
    Sr = size(trans_t,1);
    
    for k = 1:Sr
        scale = sum(trans_t(k,:));
        if scale == 0
            scale = 1;
        end
        trans_t(k,:) = trans_t(k,:)./repmat(scale,1,Sr);
    end
    tempHmm.A = trans_t;
    
    % emit
    for k = 1:Sr
        tempHmm.emit{k}.type     = 'gmm';  % a GMM with one component
        tempHmm.emit{k}.nin      = dim;
        tempHmm.emit{k}.ncentres = 1;
        tempHmm.emit{k}.priors   = 1;
        
        %%%% update covariances and mean        
        v = h3m_r.hmm{1, j}.emit{1, k}.v;
        W = h3m_r.hmm{1, j}.emit{1, k}.W;
        
        tempHmm.emit{k}.centres = h3m_r.hmm{1, j}.emit{1, k}.m;

        switch(covmode)
            case 'diag'
                if (v > dim+1)
                    tC = (1./W)/(v-dim-1);
                    Cov_new = (tC +tC)/2;
                else
                    tC = (1./W)/(v);
                    Cov_new = (tC +tC)/2;
                end
                
                tempHmm.emit{k}.covars = Cov_new;
                tempHmm.emit{k}.covar_type = 'diag';
                
            case 'full'
                if (v > dim+1)
                    tC = inv(W)/(v-dim-1);
                    Cov_new = (tC +tC')/2;
                else
                    tC = inv(W)/(v);
                    Cov_new = (tC +tC')/2;
                end
                tempHmm.emit{k}.covars = Cov_new;                
                tempHmm.emit{k}.covar_type = 'full';
                
            otherwise
                error('NOT SURPORT')
        end
 
    end
    
    h3m_new.hmm{j} = tempHmm;
 end
%  
omega_s = sum(h3m_r.alpha);
omega = h3m_r.alpha ./ omega_s;
h3m_new.omega = omega;
