function [LB, d_LB,Lsave] = vbhemh3m_colb(h3m_stats, h3m_rs, vbhemopt, bfopt)
% (internal function)
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-07-29
% Antoni B. Chan, Janet H. Hsiao, Lan Hui
% City University of Hong Kong, University of Hong Kong

% 2021-07-29: v0.80 - Hui - initial version

Nstimuli = size(h3m_rs,1);
%Z        = h3m_stats.Z;
NL       = h3m_stats.NL;
hat_Z    = h3m_stats.hat_Z;
Sum_s_L_elbo   = h3m_stats.Sum_s_L_elbo;
logOmegaTilde = h3m_stats.logOmegaTilde;
Srs      = h3m_stats.Srs;

dim      = length(vbhemopt.m0);
Kb = size(hat_Z,1);
Kr = size(hat_Z,2);

% load the infomations
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
W0mode   = h3m_rs(1).W0mode;


switch(W0mode)
    case 'iid'
        logdetW0inv = dim*log(W0inv(1,1));
    case 'diag'
        logdetW0inv = sum(log(diag(W0inv)));
    otherwise
        error('bad W0mode');
end


%Lt1 = zeros(Nstimuli,1);
% Lt2 = zeros(Nstimuli,1);
%  Lt3 = zeros(Nstimuli,1);
Lt4 = zeros(Nstimuli,1);
Lt5 = zeros(Nstimuli,1);
Lt6 = zeros(Nstimuli,1);
% Lt7 = zeros(Nstimuli,1);
%  Lt8 = zeros(Nstimuli,1);
Lt9 = zeros(Nstimuli,1);
Lt10 = zeros(Nstimuli,1);


Lt1 = sum(sum(hat_Z.*Sum_s_L_elbo));

% E[log p(Z|Omega)]  term 2
Lt2 = NL*logOmegaTilde';

logCalpha0 = gammaln(Kr*alpha0) - Kr*gammaln(alpha0);
% E[log p(Omega)]   Bishop (10.73)    term 3
Lt3 = (logCalpha0 + (alpha0-1)*sum(logOmegaTilde));  %Nstimuli*

%E [log q(z)]   term 7
Lt7 = sum(sum(hat_Z.*log(hat_Z)));
%E [log q(Omega)]  term 8
alpha = h3m_rs(1).alpha;
logCalpha = gammaln(sum(alpha)) - sum(gammaln(alpha));
Lt8 = logCalpha + (alpha -1)*logOmegaTilde';

AlogLambdaTilde = cell(1,Nstimuli);
AlogATilde      = cell(1,Nstimuli);
AlogPiTilde     = cell(1,Nstimuli);
AmWm            = cell(1,Nstimuli);
    
for s =1:Nstimuli
        
    h3m_r = h3m_rs(s);
    Kr = h3m_r.K;
    Sr = Srs(s);

    logLambdaTilde = zeros(Sr,Kr);
    logATilde      = zeros(Sr,Sr,Kr);
    logPiTilde     = zeros(Sr,Kr);   
      
    if nargin<4
     
        for j = 1: Kr
            for k = 1:Sr
            logLambdaTilde(k,j) = h3m_r.hmm{1, j}.emit{1, k}.logLambdaTilde;
            end
            logATilde(:,:,j) = h3m_r.hmm{j}.logATilde;
            logPiTilde(:,j) = h3m_r.hmm{j}.logPiTilde;
        end
        
        AlogLambdaTilde{s} = logLambdaTilde;
        AlogATilde{s}      = logATilde;
        AlogPiTilde{s}     = logPiTilde;
              
    else
        [Kb, Kr] = size(bfopt);
        for j = 1: Kr
            logLambdaTilde(:,j) = bfopt{1,j}.logLambdaTilde;
            logATilde(:,:,j) = bfopt{1,j}.logATilde;
            logPiTilde(:,j) = bfopt{1,j}.logPiTilde;  
        end       
    end   
    
     %% calculate constants
    %constants
    logCeta0 = gammaln(Sr*eta0) - Sr*gammaln(eta0);
  
    logCepsilon0 = repmat((gammaln(Sr*epsilon0) - Sr*gammaln(epsilon0)),[1,Sr]);
  
    logB0 = (v0/2)*logdetW0inv - (v0*dim/2)*log(2) ...
          - (dim*(dim-1)/4)*log(pi) - sum(gammaln(0.5*(v0+1 -[1:dim])));

    const2 = dim*log(lambda0/(2*pi));

   %% calculate each term in the LB

    
    % E[log p(pi)]   Bishop (10.73)    term 6
    Lt4(s) = Kr*logCeta0 + (eta0-1)*sum(sum(logPiTilde));

    % E[log p(A)] = sum sum E[log p(a^l_j)]   (equivalent to Bishop 10.73)
    % term 4 
    Lt5(s) = Kr*sum(logCepsilon0) + (epsilon0 -1)*sum(sum(sum(logATilde)));
  
    H = zeros(Kr,1);
    mWm = zeros(Kr,Sr);
    trW0invW = zeros(Kr,Sr);
    for j= 1:Kr
        
        lambda = zeros(Sr,1);
        v = zeros(Sr,1);
        W = zeros(dim,dim,Sr);
        m =zeros(dim,Sr);
        
        for n = 1:Sr
        lambda(n) = h3m_r.hmm{1, j}.emit{1, n}.lambda;
        v(n) =h3m_r.hmm{1, j}.emit{1, n}.v;
        W(:,:,n) = h3m_r.hmm{1, j}.emit{1, n}.W; 
        m(:,n) = h3m_r.hmm{1, j}.emit{1, n}.m'; 
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
        
        AmWm{s} = mWm;
               
         % E[log p(mu, Lambda)] term 6
        Lt61 = 0.5*sum(const2 + logLambdaTilde(:,j) - dim*lambda0./lambda - lambda0*v.*mWm(j,:)');
        Lt62 = Sr*logB0 + 0.5*(v0-dim-1)*sum(logLambdaTilde(:,j)) - 0.5*sum(v.*trW0invW(j,:)');
        Lt6(s) = Lt6(s)+Lt61+Lt62; 

        % E[log q(pi^l)] term 9
        Lt9a = logCeta + (eta-1)'*logPiTilde(:,j);

        % E[log q(a^l_j)] term 9B
        Lt9b = logCepsilon + sum((epsilon - 1).*logATilde(:,:,j),2)';
        
        Lt9(s) = Lt9(s) + Lt9a + sum(Lt9b);

        % E[ log q(mu,A)] term 10
        Lt10(s) = Lt10(s) + 0.5*sum(logLambdaTilde(:,j) + dim*log(lambda/(2*pi)) ) - 0.5*dim*Sr - H(j); 

    end
    
end
    
%% Lower-bound value
% sum all terms together
    LB = sum(Lt1) + Lt2  + Lt3 + sum(Lt4) + sum(Lt5) + sum(Lt6) - Lt7 - Lt8 - sum(Lt9) - sum(Lt10)  ;
 %   L_try = Lt1_try + Lt2  + Lt3 + Lt4 + Lt5 + Lt6 - Lt7 - Lt8 - Lt9 - Lt10 ;
 
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
     
    
%% Calculate the derivatives %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (do_deriv)     
  %% alpha0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  d_logC0.d_alpha0 = Kr*psi(0,Kr*alpha0) - Kr*psi(0,alpha0);

  sumlogOmegaTilde = sum(logOmegaTilde);

  d_Lt.d_alpha0    = Nstimuli*(d_logC0.d_alpha0 + sumlogOmegaTilde);
  
  %% eta0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  sumlogPiTilde = 0;
  d_logC0.d_eta0 = 0;
  for s = 1:Nstimuli  
      Sr = Srs(s); 
      logPiTilde = AlogPiTilde{s};
      
      d_logC0.d_eta0 = d_logC0.d_eta0 + Sr*psi(0,Sr*eta0) - Sr*psi(0,eta0); % Sr may be different    
      sumlogPiTilde  = sumlogPiTilde + sum(logPiTilde(:));
  end

  d_Lt.d_eta0  = Kr * d_logC0.d_eta0 + sumlogPiTilde;
  
  %% epsilon0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  sumlogATilde = 0;
  d_logC0.d_epsilon0 = 0;
  for s = 1: Nstimuli      
      Sr = Srs(s);
      logATilde = AlogATilde{s};
      
      d_logC0.d_epsilon0 = d_logC0.d_epsilon0 +Sr*(Sr*psi(0,Sr*epsilon0) - Sr*psi(0,epsilon0));     
      sumlogATilde = sumlogATilde + sum(logATilde(:));
  end

  d_Lt.d_epsilon0 = Kr*d_logC0.d_epsilon0 + sumlogATilde;
  
  %% v0 (nu0) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  d_logB0.d_v0  = 0.5*logdetW0inv - (dim/2)*log(2) - 0.5*sum(psi(0, 0.5*(v0+1-[1:dim])));
  tmp = 0;
  sumlogLambdaTilde = 0;
  for s = 1: Nstimuli
      Sr = Srs(s);
      logLambdaTilde = AlogLambdaTilde{s};
      
      tmp = tmp + Sr*d_logB0.d_v0;
      sumlogLambdaTilde = sumlogLambdaTilde + sum(logLambdaTilde(:));
  end
  
   d_Lt.d_v0     = Kr*tmp + 0.5*sumlogLambdaTilde;
  
  %% lambda0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  d_Lt.d_lambda0 = 0;
  for s = 1: Nstimuli
      h3m_r = h3m_rs(s);
      mWm = AmWm{s};
      Sr = Srs(s);
      for j = 1: Kr
          
          lambdaj = zeros(Sr,1);
          vj = zeros(Sr,1);
          for n = 1:Sr
          lambdaj(n,1) = h3m_r.hmm{1, j}.emit{1, n}.lambda;
          vj(n,1) = h3m_r.hmm{1, j}.emit{1, n}.v;
          end
          d_Lt.d_lambda0 = d_Lt.d_lambda0 + 0.5*sum(dim/lambda0 - dim./lambdaj(:) - vj(:).*mWm(j,:)');
      end
  end

  %% W0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   switch(W0mode)    
       case 'iid'
           d_Lt.d_W0 = 0;
           for s = 1 : Nstimuli
               Sr = Srs(s);
               h3m_r = h3m_rs(s);
               
               myW0inv = W0inv(1,1);
               myW0    = 1/myW0inv;
               d_logB0.d_W0  = -0.5*v0 * dim*myW0inv;
               d_trW0invW.d_W0 = zeros(Sr,Kr);
               
               for j =1:Kr

                   for n = 1:Sr
                       Wj = h3m_r.hmm{1, j}.emit{1, n}.W;
                       vj = h3m_r.hmm{1, j}.emit{1, n}.v;

                       d_trW0invW.d_W0(n,j) = -vj*(myW0inv^2) * trace(Wj);
                   end
               end
               
               d_Lt.d_W0 = d_Lt.d_W0 + Kr*Sr*d_logB0.d_W0 - 0.5*sum(d_trW0invW.d_W0(:));              
           end
           
        case 'diag'
            error('have not change code');
            d_Lt.d_W0  = zeros(dim,1);
            myW0inv = diag(W0inv);
            myW0    = 1./myW0inv;
            d_logB0.d_W0  = -0.5 * v0 * myW0inv;
            d_trW0invW.d_W0 = zeros(dim,Kr,Sr);
            for j = 1:Kr

                   for n = 1:Sr
                       Wj = h3m_r.hmm{1, j}.emit{1, n}.W;
                       vj = h3m_r.hmm{1, j}.emit{1, n}.v;

                    d_trW0invW.d_W0(:,k,j) = -vj*(myW0inv.^2) .*diag(Wj);
                end
       
            end
      
            d_Lt.d_W0(:) = Kr*S*d_logB0.d_W0 - 0.5*sum(sum( d_trW0invW.d_W0,2),3);
            
        otherwise
          error('bad W0mode');
   end
   
  
  %% m0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tmp =  zeros(dim,1);
    for s = 1 : Nstimuli
        Sr = Srs(s);
        h3m_r = h3m_rs(s);     
        for j =1:Kr
            for n = 1:Sr
                Wj = h3m_r.hmm{1, j}.emit{1, n}.W;
                mj = h3m_r.hmm{1, j}.emit{1, n}.m';
                vj = h3m_r.hmm{1, j}.emit{1, n}.v;

                tmp = tmp + lambda0*vj*Wj*(mj-m0);
            end
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

   