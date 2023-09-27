function [LB, d_LB,Lsave] = vbhemh3m_lb(h3m_stats, h3m_r, vbhemopt)
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-07-29
% Antoni B. Chan, Janet H. Hsiao, Lan Hui
% City University of Hong Kong, University of Hong Kong

% 2021-07-29: v0.80 - initial version

Kr = h3m_r.K;
dim      = length(vbhemopt.m0);
Sr       = vbhemopt.S;

Z        = h3m_stats.Z;
Nj       = h3m_stats.Nj;
hat_Z    = h3m_stats.hat_Z;
L_elbo   = h3m_stats.L_elbo ;
logOmegaTilde = h3m_stats.logOmegaTilde;



logLambdaTilde = zeros(Sr,Kr);
logATilde      = zeros(Sr,Sr,Kr);
logPiTilde     = zeros(Sr,Kr);

for j = 1: Kr
    for n= 1:Sr
        logLambdaTilde(n,j) =h3m_r.hmm{1, j}.emit{1, n}.logLambdaTilde;
    end
    logATilde(:,:,j) = h3m_r.hmm{j}.logATilde;
    logPiTilde(:,j) = h3m_r.hmm{j}.logPiTilde;
end



if isfield(h3m_stats, 'do_deriv')
    do_deriv = h3m_stats.do_deriv;
else
    do_deriv = 0;
end

clipped  = vbhemopt.hyp_clipped;

alpha0   = vbhemopt.alpha0;
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
W0mode   = h3m_r.W0mode;




%% calculate constants
switch(W0mode)
    case 'iid'
        logdetW0inv = dim*log(W0inv(1,1));
    case 'diag'
        logdetW0inv = sum(log(diag(W0inv)));
    otherwise
        error('bad W0mode');
end

%constants
logCalpha0 = gammaln(Kr*alpha0) - Kr*gammaln(alpha0);
logCeta0 = gammaln(Sr*eta0) - Sr*gammaln(eta0);

logCepsilon0 = repmat((gammaln(Sr*epsilon0) - Sr*gammaln(epsilon0)),[1,Sr]);

logB0 = (v0/2)*logdetW0inv - (v0*dim/2)*log(2) ...
    - (dim*(dim-1)/4)*log(pi) - sum(gammaln(0.5*(v0+1 -[1:dim])));

const2 = dim*log(lambda0/(2*pi));

alpha = h3m_r.alpha;
logCalpha = gammaln(sum(alpha)) - sum(gammaln(alpha));

%% calculate each term in the LB

Lt1 = sum(sum(Z.*L_elbo));

% E[log p(Z|Omega)]  term 2
Lt2 = Nj*logOmegaTilde';

% E[log p(pi)]   Bishop (10.73)    term 3
Lt3 = Kr*logCeta0 + (eta0-1)*sum(sum(logPiTilde));

% E[log p(A)] = sum sum E[log p(a^l_j)]   (equivalent to Bishop 10.73)
% term 4
Lt4 = Kr*sum(logCepsilon0) + (epsilon0 -1)*sum(sum(sum(logATilde)));

Lt5 = 0;
% E[log p(Omega)]   Bishop (10.73)    term 6
Lt6 = logCalpha0 + (alpha0-1)*sum(logOmegaTilde);

%E [log q(z)]   term 7
Lt7 = sum(sum(hat_Z.*log(hat_Z)));

%E [log q(Omega)]  term 8
Lt8 = logCalpha + (alpha -1)*logOmegaTilde';

Lt9 = 0;
Lt10 = 0;
%   Lt11 = sum(sum(Z.*L_x));

H = zeros(Kr,1);
mWm = zeros(Kr,Sr);
trW0invW = zeros(Kr,Sr);
for j= 1:Kr
    
    for k = 1:Sr
        lambda(k,1) = h3m_r.hmm{1, j}.emit{1, k}.lambda;
        v(k,1) =  h3m_r.hmm{1, j}.emit{1, k}.v;
        m(:,k) = h3m_r.hmm{1, j}.emit{1, k}.m;
        switch(vbhemopt.emit.covar_type)
            case 'diag'
                W(:,:,k) = diag(h3m_r.hmm{1, j}.emit{1, k}.W);
            case 'full'
                W(:,:,k) = h3m_r.hmm{1, j}.emit{1, k}.W;
        end
    end
    
    eta = h3m_r.hmm{j}.eta;
    epsilon = h3m_r.hmm{j}.epsilon;
    
    logCeta = gammaln(sum(eta)) - sum(gammaln(eta));
    logCepsilon = zeros(1,Sr);
    
    for k = 1:Sr
        logCepsilon(k) = gammaln(sum(epsilon(k,:))) - sum(gammaln(epsilon(k,:)));
    end
    
    
    for k = 1:Sr
        logBk = -(v(k)/2)*log(det(W(:,:,k))) - (v(k)*dim/2)*log(2)...
            - (dim*(dim-1)/4)*log(pi) - sum(gammaln(0.5*(v(k) + 1 - [1:dim])));
        H(j) = H(j) - logBk - 0.5*(v(k) - dim - 1)*logLambdaTilde(k,j) + 0.5*v(k)*dim;
        diff = m(:,k) - m0;
        mWm(j,k) = diff'*W(:,:,k)*diff;
        trW0invW(j,k) = trace(W0inv*W(:,:,k));
    end

    % E[log p(mu, Lambda)] term 5
    Lt51 = 0.5*sum(const2 + logLambdaTilde(:,j) - dim*lambda0./lambda - lambda0*v.*mWm(j,:)');
    Lt52 = Sr*logB0 + 0.5*(v0-dim-1)*sum(logLambdaTilde(:,j)) - 0.5*sum(v.*trW0invW(j,:)');
    Lt5 = Lt5+Lt51+Lt52;
    
    % E[log q(pi^l)] term 9
    Lt9a = logCeta + (eta-1)'*logPiTilde(:,j);
    
    % E[log q(a^l_j)] term 9B
    Lt9b = logCepsilon + sum((epsilon - 1).*logATilde(:,:,j),2)';
    
    if 0
        for k = 1:Sr
            xLt9b(k) = logCepsilon(k) + (epsilon(k,:) - 1)*logATilde(k,:,j)';
        end
        err = [];
        err(end+1) = sum(abs(xLt9b - Lt9b));
        
        if any(err>1e-10)
            warning('mismatch error');
            keyboard
        end
    end
    
    Lt9 = Lt9 + Lt9a + sum(Lt9b);
    
    % E[ log q(mu,A)] term 10
    Lt10 = Lt10 +0.5*sum( logLambdaTilde(:,j) + dim*log(lambda/(2*pi))) - 0.5*dim*Sr - H(j);
    
end

%% Lower-bound value
% sum all terms together
LB = Lt1 + Lt2  + Lt3 + Lt4 + Lt5 + Lt6 - Lt7 - Lt8 - Lt9 - Lt10  ;
%   L_try = Lt1_try + Lt2  + Lt3 + Lt4 + Lt5 + Lt6 - Lt7 - Lt8 - Lt9 - Lt10 ;
if nargout>2
Lsave.Lt1 =  Lt1;
Lsave.Lt2 =  Lt2;
Lsave.Lt3 =  Lt3;
Lsave.Lt4 =  Lt4;
Lsave.Lt5 =  Lt5;
Lsave.Lt6 =  Lt6;
Lsave.Lt7 =  Lt7;
Lsave.Lt8 =  Lt8;
Lsave.Lt9 =  Lt9;
Lsave.Lt10 =  Lt10;
end

%% Calculate the derivatives %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (do_deriv)
    %% alpha0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    d_logC0.d_alpha0 = Kr*psi(0,Kr*alpha0) - Kr*psi(0,alpha0);
    
    sumlogOmegaTilde = sum(logOmegaTilde);
    
    d_Lt.d_alpha0    = d_logC0.d_alpha0 + sumlogOmegaTilde;
    
    %% eta0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    d_logC0.d_eta0 = Sr*psi(0,Sr*eta0) - Sr*psi(0,eta0);
    
    sumlogPiTilde  = sum(logPiTilde(:));
    
    d_Lt.d_eta0  = Kr *d_logC0.d_eta0 + sumlogPiTilde;
    
    %% epsilon0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    d_logC0.d_epsilon0 = Sr*psi(0,Sr*epsilon0) - Sr*psi(0,epsilon0);
    
    sumlogATilde = sum(logATilde(:));
    
    d_Lt.d_epsilon0    = Kr*Sr*d_logC0.d_epsilon0 + sumlogATilde;
    
    %% v0 (nu0) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    d_logB0.d_v0  = 0.5*logdetW0inv - (dim/2)*log(2) - 0.5*sum(psi(0, 0.5*(v0+1-[1:dim])));
    d_Lt.d_v0     = Kr*Sr*d_logB0.d_v0 + 0.5*sum(logLambdaTilde(:));
    
    %% lambda0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    d_Lt.d_lambda0 = 0;
    
    for j = 1: Kr       
        for k = 1:Sr
            lambda(k) = h3m_r.hmm{1, j}.emit{1, k}.lambda;
            v(k) =  h3m_r.hmm{1, j}.emit{1, k}.v;            
        end        
        d_Lt.d_lambda0 = d_Lt.d_lambda0 + 0.5*sum(dim/lambda0 - dim./lambda(:) - v(:).*mWm(j,:)');
    end
        
    %% W0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch(W0mode)
        case 'iid'
            
            myW0inv = W0inv(1,1);
            myW0    = 1/myW0inv;
            d_logB0.d_W0  = - 0.5*v0*dim*myW0inv;
            
            d_trW0invW.d_W0 = zeros(Sr,Kr);
            for j = 1:Kr
                for k= 1:Sr
                    v =  h3m_r.hmm{1, j}.emit{1, k}.v;
                    W= h3m_r.hmm{1, j}.emit{1, k}.W;

                    switch(vbhemopt.emit.covar_type)
                        case 'diag'
                            d_trW0invW.d_W0(k,j) = -v*(myW0inv^2) * sum(W);
                        case 'full'
                            d_trW0invW.d_W0(k,j) = -v*(myW0inv^2) * trace(W);                          
                    end                   
                end
            end
            
            d_Lt.d_W0 = Kr*Sr*d_logB0.d_W0 - 0.5*sum(d_trW0invW.d_W0(:));
            
            
        case 'diag'
            
            d_Lt.d_W0  = zeros(dim,1);
            myW0inv = diag(W0inv);
            myW0    = 1./myW0inv;
            d_logB0.d_W0  = -0.5 * v0 * myW0inv;
            d_trW0invW.d_W0 = zeros(dim,Kr,Sr);
            for j = 1:Kr
                for k=1:Sr
                    v =  h3m_r.hmm{1, j}.emit{1, k}.v;
                    W= h3m_r.hmm{1, j}.emit{1, k}.W;
                    switch(vbhemopt.emit.covar_type)
                        case 'diag'
                            d_trW0invW.d_W0(:,k,j) = -v*(myW0inv.^2) .*(W)';
                        case 'full'
                            d_trW0invW.d_W0(:,k,j) = -v*(myW0inv.^2) .*diag(W);
                    end
                    
                end
                
            end
            
            d_Lt.d_W0(:) = Kr*Sr*d_logB0.d_W0 - 0.5*sum(sum( d_trW0invW.d_W0,2),3);
            
        otherwise
            error('bad W0mode');
    end
    
    
    %% m0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tmp =  zeros(dim,1);
    for j =1:Kr        
        for k=1:Sr          
            v =  h3m_r.hmm{1, j}.emit{1, k}.v;
            m = h3m_r.hmm{1, j}.emit{1, k}.m;
            W = h3m_r.hmm{1, j}.emit{1, k}.W;
            tmp = tmp + lambda0*v*W*(m'-m0);
        end
    end
    d_Lt.d_m0 = tmp;
    
    
    %% set derivatives to 0 when moving beyond extreme values of hyps will increase LL:
    % 1) if at max hyp, prevent increasing hyp if LL will increase
    % 2) if at min hyp, prevent decreasing hyp if LL will increase
    
    hnames = fieldnames(clipped);
    for i=1:length(hnames)
        myhname = hnames{i};
        mydname = ['d_' myhname];
        
        for j=1:length(clipped.(myhname))
            if (clipped.(myhname)(j) == +1) && (d_Lt.(mydname)(j) > 0)
                d_Lt.(mydname)(j) = 0;
            end
            if (clipped.(myhname)(j) == -1) && (d_Lt.(mydname)(j) < 0)
                d_Lt.(mydname)(j) = 0;
            end
        end
    end
    
    %% calculate derivatives of log (transformed) hyps
    d_LB.d_logalpha0   = d_Lt.d_alpha0   * alpha0;
    d_LB.d_logeta0     = d_Lt.d_eta0   * eta0;
    d_LB.d_logepsilon0 = d_Lt.d_epsilon0 * epsilon0;
    d_LB.d_logv0D1     = d_Lt.d_v0 * (v0-dim+1);
    d_LB.d_sqrtv0D1    = d_Lt.d_v0 * 2 * sqrt(v0-dim+1);
    d_LB.d_loglambda0    = d_Lt.d_lambda0 * lambda0;
    d_LB.d_sqrtlambda0   = d_Lt.d_lambda0 * 2 * sqrt(lambda0);
    d_LB.d_sqrtW0inv   = bsxfun(@times, d_Lt.d_W0, (myW0.^1.5) * (-2));
    d_LB.d_logW0       = bsxfun(@times, d_Lt.d_W0, myW0);
    d_LB.d_m0          = d_Lt.d_m0;
    
    % debugging
    %d_LB.LBterms = [Lt1 Lt2 Lt3 Lt4 Lt5 Lt6 Lt7 Lt8];
    
else
    d_LB = [];
    
end

   