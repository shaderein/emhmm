function [hmm_new,Syn_stats]= vbhem_mstep_component(estats, vbhemopt)
% VBHEM-MSTEP
% update the component/hmm of h3m_r
% update hyperparameters not parameters
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-07-29
% Antoni B. Chan, Janet H. Hsiao, Lan Hui
% City University of Hong Kong, University of Hong Kong

% 2021-07-29: v0.80 - initial version

Sr  = estats.Sr;
dim = estats.dim;
covmode = vbhemopt.emit.covar_type;

eta0     = vbhemopt.eta0;
epsilon0 = vbhemopt.epsilon0;
m0       = vbhemopt.m0;
lambda0  = vbhemopt.lambda0;

if numel(vbhemopt.W0) == 1       
    % isotropic W
    W0  = vbhemopt.W0*eye(dim);  
else
 % diagonal W
    if numel(vbhemopt.W0) ~= dim
        error(sprintf('vbhemopt.W should have dimension D=%d for diagonal matrix', dim));
    end

    W0  = diag(vbhemopt.W0);
end

if vbhemopt.v0<=dim-1       
    error('v0 not large enough');
end
v0       = vbhemopt.v0; % should be larger than p-1 degrees of freedom (or dimensions)
W0inv    = inv(W0);
     


[Nj_rho1,Nj_rho2rho,Nj_rho,y_bar,S_plus_C]= vbhem_compute_Statistics(estats, vbhemopt);
Syn_stats.Nj_rho1 = Nj_rho1;
Syn_stats.Nj_rho2rho = Nj_rho2rho;
Syn_stats.Nj_rho = Nj_rho;
Syn_stats.y_bar = y_bar;
Syn_stats.S_plus_C = S_plus_C;

hmm_new.eta = eta0 + Nj_rho1;
hmm_new.epsilon = epsilon0 + Nj_rho2rho;

hmm_new.emit = {};

for k = 1:Sr
    
    hmm_new.emit{k}.lambda = lambda0 + Nj_rho(k);
    hmm_new.emit{k}.v = v0 + Nj_rho(k) +1;
    
    m_new = (lambda0*m0 + Nj_rho(k)*y_bar(k,:)')/(lambda0 + Nj_rho(k));
    hmm_new.emit{k}.m =  m_new';
    
    mult1 = lambda0*Nj_rho(k)/(lambda0 + Nj_rho(k));
    diff3 = y_bar(k,:) - m0';
    
    switch covmode
        case 'diag'
            tW = inv( W0inv + Nj_rho(k)*(diag(S_plus_C(k,:)))...
                + mult1*(diff3'*diff3) );
            W_new = (tW+ tW')/2;
            hmm_new.emit{k}.W =  diag(W_new)';
            
        case 'full'
            tW = inv( W0inv + Nj_rho(k)* reshape(S_plus_C(k,:,:),[dim,dim])...
                + mult1*(diff3'*diff3) );
            W_new = (tW+ tW')/2;
            hmm_new.emit{k}.W =  W_new;
    end


end



