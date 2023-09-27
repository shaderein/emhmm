function [Nj_rho1,Nj_rho2rho,Nj_rho,y_bar,S_plus_C]= vbhem_compute_Statistics(estats, vbhmopt)
% compute the sythetic statistic in the vbhem
% Nj for w/ alpha
% Nj_rho1 for pi/eat
% Nj_rho2rho for A/epsilon
% Nj_rho , y_bar_jrho, S_jrho, C_jrho for mean & cov 
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-07-29
% Antoni B. Chan, Janet H. Hsiao, Lan Hui
% City University of Hong Kong, University of Hong Kong

% 2021-07-29: v0.80 - initial version


Sr  = estats.Sr;
dim = estats.dim;
Kb  = estats.Kb;
covmode = vbhmopt.emit.covar_type;

Z_Ni = estats.Z_Ni;


%% compute the summary statistics
Nj_rho1   = zeros(Sr,1);
Nj_rho2rho   = zeros(Sr,Sr);
Nj_rho     = zeros(Sr,1);
y_bar    = zeros(Sr,dim);

switch covmode
    case 'diag'
        S_plus_C     = zeros(Sr,dim);
        
    case 'full'
        S_plus_C     = zeros(Sr,dim,dim);
end



% loop all the components of the base mixture
for i = 1 : Kb
    
    if (Z_Ni(i)>1e-8)
        
        nu     = estats.nu_1{i};           % this is a 1 by N vector
        xi     = estats.sum_xi{i};         % this is a N by N matrix (from - to)
        up_pr  = estats.emit_pr{i}; % this is a N by M matrix
        up_mu  = estats.emit_mu{i}; % this is a N by dim by M matrix
        up_M  = estats.emit_M{i};  % this is a N by dim by M matrix [diagonal covariance]
        
        
        Nj_rho1   = Nj_rho1 + Z_Ni(i) * nu';      % update eta
        Nj_rho2rho = Nj_rho2rho + Z_Ni(i) * xi;     % update epsilon
        Nj_rho     = Nj_rho + Z_Ni(i) * up_pr;
        
        y_bar = y_bar + Z_Ni(i) * up_mu;
        
        S_plus_C = S_plus_C + Z_Ni(i) * up_M;

    else
        continue
    end
end

Nj_rho    = Nj_rho + 1e-50;

y_bar = y_bar./(Nj_rho*ones(1,dim));

switch covmode
    case 'diag'
        for k=1:Sr
            S_plus_C(k,:) = S_plus_C(k,:)./(Nj_rho(k)*ones(1,dim));

            S_plus_C(k,:) =  S_plus_C(k,:) - y_bar(k,:).*y_bar(k,:);
        end
        
    case 'full'
        
        
        for k=1:Sr
            S_plus_C(k,:,:) = S_plus_C(k,:,:)./(Nj_rho(k)*ones(1,dim));
        
            S_plus_C(k,:,:) =  reshape(S_plus_C(k,:,:),[dim,dim]) - y_bar(k,:)'*y_bar(k,:);
        end
        
end

if (Sr==1)
    Nj_rho2rho  = 1e-12;
end



