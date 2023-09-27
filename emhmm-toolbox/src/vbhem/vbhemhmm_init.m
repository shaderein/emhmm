function [h3m_r,vbhemopt] = vbhemhmm_init(h3m_b,vbhemopt)
% vbhemhmm_init - initialize vbh3m (internal function)
% the idea for initialization the hyperparameters is based on the
% initialization ways in the VHEM, using the initial h3m_r to update the count number
% count number has N^j, N^j_rho1,N^j_(rho,rho'), N^j_rho'
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-07-29
% Antoni B. Chan, Janet H. Hsiao, Lan Hui
% City University of Hong Kong, University of Hong Kong

% 2020-9-23: the new version the m in the emit is [m_1,m_2]; the old type
% is[m1;m2]
% 2021-07-29: v0.80 - initial version

VERBOSE_MODE = vbhemopt.verbose;

h3m   = {};
Kb    = h3m_b.K;
Kr    = vbhemopt.K;
h3m.K = Kr;
Sr     = vbhemopt.S;   
N     = vbhemopt.S;
h3m.S = Sr;
Nv    = vbhemopt.Nv *Kb;
T     = vbhemopt.tau;
dim   = h3m_b.hmm{1}.emit{1}.nin;

[vbhemopt, hyp_clipped] = vbhem_clip_hyps(vbhemopt);
vbhemopt.hyp_clipped = hyp_clipped;
    
% setup hyperparameters
alpha0   = vbhemopt.alpha0;
eta0     = vbhemopt.eta0;
epsilon0 = vbhemopt.epsilon0;
m0       = vbhemopt.m0;
lambda0  = vbhemopt.lambda0;

if numel(vbhemopt.W0) == 1
  % isotropic W
  W0       = vbhemopt.W0*eye(dim);  
  h3m.W0mode   = 'iid';

else
  % diagonal W
  if numel(vbhemopt.W0) ~= dim
    error(sprintf('vbhemopt.W should have dimension D=%d for diagonal matrix', dim));
  end
  W0       = diag(vbhemopt.W0);
  h3m.W0mode   = 'diag';
end


if vbhemopt.v0<=dim-1
  error('v0 not large enough');
end
v0       = vbhemopt.v0; % should be larger than p-1 degrees of freedom (or dimensions)
W0inv    = inv(W0);


switch vbhemopt.initmode

    % Nv hase times Kb
    case 'baseem'
        
        NLr= (Nv/Kr);
        
        for j=1:Kr
            
            for n=1:Sr    % state
                
                if isfield(vbhemopt, 'prerand')
                    
                    randomb = vbhemopt.prerand{1,1}(j,n);
                    randomg = vbhemopt.prerand{2,1}(j,n);
                    
                else
                    randomb = randi(Kb);
                    randomg = randi(length(h3m_b.hmm{randomb}.emit));
                end
                
                h3m.hmm{j}.emit{n}.lambda = lambda0 + (NLr)/Sr ;
                h3m.hmm{j}.emit{n}.v = v0 + (NLr)/Sr+1;
                h3m.hmm{j}.emit{n}.m= h3m_b.hmm{randomb}.emit{randomg}.centres;
                                        
                switch(h3m_b.hmm{1}.emit{1}.covar_type)
                    case 'diag'
                        tmp =  1./((h3m.hmm{j}.emit{n}.v-dim-1) * (h3m_b.hmm{randomb}.emit{randomg}.covars));
                        h3m.hmm{j}.emit{n}.W = tmp;
                        
                    case 'full'
                        h3m.hmm{j}.emit{n}.W= inv((h3m.hmm{j}.emit{n}.v-dim-1) * (h3m_b.hmm{randomb}.emit{randomg}.covars));     
                end

            end
            [prior, A] = makeAprior(Sr, vbhemopt.initopt.mode);
            h3m.hmm{j}.eta   = prior * NLr + eta0 ;  %initial uniform P(pi_rho1) = 1/Kr
            h3m.hmm{j}.epsilon = (A * NLr)/(Sr) +epsilon0;
            
        end
        
        omega = rand(Kr,1);
        omega = omega/sum(omega);
        h3m.LogL = -inf;
        h3m.alpha = alpha0 + omega'*Nv;
        h3m_r =h3m;

    % use hem-g3m    
    case {'gmmNew', 'gmmNew2'}
        
    
        % get all the emission distributions, and using the VHEM to reduce
        % them to a GMM with N (gmmNew) or N*Kr components.
        % if gmmNew - estimate N Gaussians for N states - share Gaussians across groups
        % if gmmNew2 - estimate N*Kr Gaussians - each group has different Gaussians.
        % the transition and prior is random for each hmm. (except the 'm'
        % case)
    
        virtualSamples = Nv ;
        iterations = vbhemopt.initopt.iter;
        % fit a gaussian mixture to the data...
    
        % 2017-05-30: ABC
        %% NEW CODE: handles different number of ROIs
        getN = @(x) length(x.prior);
        numbc = sum(cellfun(getN, h3m_b.hmm));

        gmms_all = cell(1,numbc);
        gmms_hmm_index = zeros(1,numbc);
        gmms_hmm_emit_index = zeros(1,numbc);
    
        curj = 1;
        for i = 1 : h3m_b.K
          myS = length(h3m_b.hmm{i}.prior);
          newind = curj:(curj+myS-1);
          gmms_all( curj:(curj+myS-1) ) = h3m_b.hmm{i}.emit;

          gmms_hmm_index(newind) = i;
          gmms_hmm_emit_index(newind) = 1:length(newind);

          curj = curj+myS;
        end

    
        trainMix(length(gmms_all)) = struct(gmms_all{1});
        for pippo = 1 : length(gmms_all)
          trainMix(pippo).mix = gmms_all{pippo};

        end
    
        if strcmp(vbhemopt.initmode, 'gmmNew2')
          tmpK = N*Kr; % separate ROIs for each group
          useGauss = reshape(randperm(tmpK),Kr,N);  % [Kr x N] index into ROIs
        else
          tmpK = N; % share ROIs for each group
          useGauss = repmat(1:N,Kr,1);              % [Kr x N] index into ROIs
        end
    
        % ABC 2017-05-18: BUG fix: should be N (usually Kr>N for tagging)
        [reduced_out, ~, lp_out] = GMM_MixHierEM(trainMix, tmpK, virtualSamples, iterations);
        %[reduced_out] = GMM_MixHierEM(trainMix, Kr, virtualSamples, iterations);
        % post_out is the position,
        %lp_out is like the hat_z_ij
    
    
        % map from HMM emission to group ROI
        emit_map_ll = cell(1,Kb);
        for i=1:numbc
           emit_map_ll{gmms_hmm_index(i)}(:,gmms_hmm_emit_index(i)) = lp_out(:,i);
        end
        
        % assign base HMMS to reduced HMM using emissions
        hmm_map = zeros(1,Kb);
        hmm_map_emit = cell(1,Kb);
        for i=1:Kb
          % compute sum of LL between group ROIs and HMM emissions
          emit_LL = zeros(1,Kr);
          for j=1:Kr
            emit_LL(j) = sum(sum(emit_map_ll{i}(useGauss(j,:),:)));
            % the emission of jth HMM in reduced model consist of useGauss(j,:)
          end
          % pick largest sum LL as the group
          [~, myj] = max(emit_LL);
          % map the base ROIs to group ROIs.
          [~, myii] = max(emit_map_ll{i}(useGauss(myj,:),:),[], 1);   %?
          hmm_map(i) = myj;
          hmm_map_emit{i} = myii;                       
        end
    
        % build new transition matrices
        newA = cell(1,Kr);
        newprior = cell(1,Kr);
        newN = zeros(1,Kr);   
        for j=1:Kr
          newA{j} = zeros(N,N);
          newprior{j} = zeros(N,1);
          newN(j) = 0;
        end
        for i=1:Kb      
          mymap = hmm_map_emit{i}; % get the ROI map from base to reduced
          myj = hmm_map(i);      % get the corresponding group
          % for each base hidden state
          for k=1:length(mymap)
            nk = mymap(k);
            % add probability to corresponding reduced hidden state
            newprior{myj}(nk) = newprior{myj}(nk) + h3m_b.hmm{i}.prior(k);
            for l=1:length(mymap)
              nl = mymap(l);
              newA{myj}(nk,nl) = newA{myj}(nk,nl) + h3m_b.hmm{i}.A(k,l);
            end
          end
          newN(myj) = newN(myj)+1;
        end
        % normalize
        for j=1:Kr
          newprior{j} = newprior{j} / sum(newprior{j});
          newA{j} = bsxfun(@rdivide, newA{j}, sum(newA{j},2));
        end

    
        for j = 1 : Kr

          emit = cell(1,N);
          for n = 1 : N
            myn = useGauss(j, n);        
            reduced_out_use.centres = reduced_out.centres(myn,:);
            switch(h3m_b.hmm{1}.emit{1}.covar_type)
              case 'diag'
                reduced_out_use.covars = reduced_out.covars(myn,:);
              case 'full'
                reduced_out_use.covars = reduced_out.covars(:,:,myn);
            end
            reduced_out_use.ncentres = 1;
            reduced_out_use.priors = 1;
            reduced_out_use.nin = reduced_out.nin;
            reduced_out_use.covar_type = reduced_out.covar_type;
            reduced_out_use.type = reduced_out.type;

            emit{n} = reduced_out_use;
            % and randomize priors ...
            emit{n}.priors = rand(size(emit{n}.priors));
            emit{n}.priors = emit{n}.priors / sum(emit{n}.priors);
          end

          % ABC - 2017-05-18: Make it random prior/A      
          % original case was random      
          if any(vbhemopt.initopt.mode == 'm')
            h3m.hmms{j}.prior = newprior{j};
            h3m.hmms{j}.A     = newA{j};
          else
            [prior, A] = makeAprior(N, vbhemopt.initopt.mode);
            h3m.hmms{j}.prior = prior;
            h3m.hmms{j}.A     = A;
          end  

          h3m.hmms{j}.emit  = emit;
        end
        omega = rand(1,Kr);
        omega = omega/sum(omega);
        h3m.omega = omega;
        h3m.LogL = -inf;    

        %%%using the h3m to get back the hyperparameter
        Nsj = h3m.omega* virtualSamples;
        h3m_vb.K = h3m.K;
        h3m_vb.S = h3m.S;
        h3m_vb.W0mode = h3m.W0mode;
        
        h3m_vb.hmm = {};
    
        for j=1:Kr

          h3m_vb.hmm{j}.eta   = h3m.hmms{1, j}.prior*Nsj(j) + eta0 ;  %initial uniform P(pi_rho1) = 1/Kr 
          h3m_vb.hmm{j}.epsilon     = h3m.hmms{1, j}.A *Nsj(j) +epsilon0; 
          Nsj_rho = Nsj(j)*ones(1,Sr)/Sr;                    

          h3m_vb.hmm{j}.emit = {};

          for n=1:Sr
              
              h3m_vb.hmm{j}.emit{1,n}.v = v0 + Nsj_rho(n) + 1;
              h3m_vb.hmm{j}.emit{1,n}.lambda = lambda0 + Nsj_rho(n);
              
              h3m_vb.hmm{j}.emit{1,n}.m = h3m.hmms{j}.emit{n}.centres;


              switch(h3m_b.hmm{1}.emit{1}.covar_type)
              case 'diag'
                 h3m_vb.hmm{j}.emit{n}.W = 1./((h3m_vb.hmm{j}.emit{1,n}.v-dim-1) *h3m.hmms{j}.emit{n}.covars);
              case 'full'
                 h3m_vb.hmm{j}.emit{n}.W = inv((h3m_vb.hmm{j}.emit{1,n}.v-dim-1) *((h3m.hmms{j}.emit{n}.covars)));
              end    
          end      
        end
    
        h3m_vb.alpha = alpha0 + Nsj;
        h3m_r =h3m_vb;

        
    case 'wtkmeans'      %'
        
        getN = @(x) length(x.prior);
        numbc = sum(cellfun(getN, h3m_b.hmm));
    
        % combine all the gausians as a gmm
        gmms_all = cell(1,numbc);
        alpha = zeros(1,numbc);

        curj = 1;
        for i = 1 : h3m_b.K
          myS = length(h3m_b.hmm{i}.prior);
          newind = curj:(curj+myS-1);
          gmms_all( curj:(curj+myS-1) ) = h3m_b.hmm{i}.emit;
          p = makeGMMweights(h3m_b.hmm{i}.prior', h3m_b.hmm{i}.A, vbhemopt.initopt.mode);           
          alpha(newind) = p;
          curj = curj+myS;
        end
    
        alpha = alpha/sum(alpha);
        gmms.alpha = alpha;
        gmms.mix   = gmms_all;
        
        %%%%%%%%
        
        % get all the base means
        mumtx = zeros(dim, numbc);
        for j = 1: numbc
            mumtx(:,j) = gmms.mix{1, j}.centres';  %dim*n   '
        end
        
        it_max = 100;
        % at firt, use K-means seperate means to Kr clusters, get the
        % assignment information
        
        % give a start seed for the kmeans
        rng(vbhemopt.wtseed, 'twister');
        [~,init_center] = kmeans(mumtx',Kr,'Replicates',1);    %'input is n*dim, output is k*dim
        % use Kmeans results as initial for weighted kmeans
        [cluster_idx,~] = my_weighted_kmeans(Kr,it_max,mumtx,gmms.alpha,init_center');
        
        %%%%% cluster the means in one cluster to get the members of each states
        cluster_centers = cell(Kr,1);
        empty_ind = zeros(Kr,1);
        
        for i = 1: Kr
            meani = mumtx(:,cluster_idx==i);
            
            if isempty(meani)
                cluster_centers{i} = [];
                empty_ind(i) = 1;
                continue
            end
            
            tmpK = Sr;
            N_mean = size(meani,2);
         
          if (N_mean<=Sr)

              Nextra = Sr-N_mean;
              extr_mean = repmat(meani(:,1),[1,Nextra]);
              
              cluster_centers{i} = [meani,extr_mean];
          else
              
            rng(vbhemopt.wtseed+1, 'twister');
            [~,init_center1] = kmeans(meani',tmpK,'Display','off','Replicates',1); 
          
            temp = init_center1';  %   dim*S

            if size(temp,2)~=Sr
                temp1 = temp;
                extra_s = Sr-size(temp1,2);
                
                temp = [temp1,repmat(temp1(:,1),[1,extra_s])];
                
            end
            cluster_centers{i} = temp;
           
          end
        end
        
        first_unempty = find(empty_ind~=1, 1, 'first');
        if sum(empty_ind)~=0
          for i = 1:Kr       
              if empty_ind(i)
                  cluster_centers{i} = cluster_centers{first_unempty};
              end
          end
        end
        
       
        NJ = Nv *Kb/Kr;
        h3m.hmm = {};
    
        for j=1:Kr
            
          [prior, A] = makeAprior(N, vbhemopt.initopt.mode);

          h3m.hmm{j}.eta   = prior * NJ + eta0 ;  %initial uniform P(pi_rho1) = 1/Kr 
          h3m.hmm{j}.epsilon  = A * NJ + epsilon0; 
          Nsj_rho = NJ*ones(1,Sr)/Sr;                    
%           h3m.hmm{j}.lambda = lambda0 + Nsj_rho';
%           h3m.hmm{j}.v = v0 + Nsj_rho' + 1;
          h3m.hmm{j}.emit = {};
          
          cluster_center = cluster_centers{j};
          
          if size(cluster_center,2)<Sr
              keyboard
          end
          
          for n=1:Sr

              h3m.hmm{j}.emit{n}.m = cluster_center(:,n)';
              h3m.hmm{j}.emit{n}.v = v0 + Nsj_rho(n) + 1;
              h3m.hmm{j}.emit{n}.lambda = lambda0 + Nsj_rho(n);

              switch(h3m_b.hmm{1}.emit{1}.covar_type)
              case 'diag'
                 tmp =  1./((h3m.hmm{j}.emit{n}.v-dim-1) * (h3m_b.hmm{1, 1}.emit{1, 1}.covars));
                 h3m.hmm{j}.emit{n}.W = tmp;
                 
              case 'full'
                  h3m.hmm{j}.emit{n}.W= inv((h3m.hmm{j}.emit{n}.v-dim-1) * (h3m_b.hmm{1, 1}.emit{1, 1}.covars));
                  
              end    
          end      
        end
    
        h3m.alpha = alpha0 + NJ*ones(1,Kr);
        h3m_r =h3m;
        
        
    case 'inith3m'
        
        % use the h3m we have got
        h3m.alpha = vbhemopt.inithmm.alpha;
        h3m.hmm = {};
        
        
        for j = 1: Kr
            
            varp =vbhemopt.inithmm.hmm{j}.varpar;
            
            h3m.hmm{j}.eta = varp.alpha;
            h3m.hmm{j}.epsilon = varp.epsilon;

            m = varp.m;
            W = varp.W;  
            
            Sr = size(m,2);
            h3m.hmm{j}.emit = {};
            for k = 1: Sr
                h3m.hmm{j}.emit{k}.m = m(:,k)';
                h3m.hmm{j}.emit{k}.W = W(:,:,k);
                h3m.hmm{j}.emit{k}.v = varp.v(k);
                h3m.hmm{j}.emit{k}.lambda = varp.beta(k);
            end
           
        end
        h3m_r = h3m;
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
    case {'gmmNew_r', 'gmmNew2_r'}
        % using the emission map information to generate a hat_z, generate
        % the theta and then upadate the parameters using hat_z and theta.
        
    % gmmNew - estimate N Gaussians for N states - share Gaussians across groups
    % gmmNew2 - estimate N*Kr Gaussians - each group has different Gaussians.
   
        virtualSamples = Nv ;
        iterations = vbhemopt.initopt.iter;

        Sr = N;
        % fit a gaussian mixture to the data...

        % 2017-05-30: ABC
        %% NEW CODE: handles different number of ROIs
        getN = @(x) length(x.prior);
        numbc = sum(cellfun(getN, h3m_b.hmm));

        gmms_all = cell(1,numbc);
        gmms_hmm_index = zeros(1,numbc);
        gmms_hmm_emit_index = zeros(1,numbc);

        curj = 1;
        for i = 1 : h3m_b.K
          myS = length(h3m_b.hmm{i}.prior);
          newind = curj:(curj+myS-1);
          gmms_all( curj:(curj+myS-1) ) = h3m_b.hmm{i}.emit;

          gmms_hmm_index(newind) = i;
          gmms_hmm_emit_index(newind) = 1:length(newind);

          curj = curj+myS;
        end

        trainMix(length(gmms_all)) = struct(gmms_all{1});
        for pippo = 1 : length(gmms_all)
          trainMix(pippo).mix = gmms_all{pippo};

        end
        
        if strcmp(vbhemopt.initmode, 'gmmNew2_r')
          tmpK = N*Kr; % separate ROIs for each group
          useGauss = reshape(randperm(tmpK),Kr,N);  % [Kr x N] index into ROIs
        else
          tmpK = N; % share ROIs for each group
          useGauss = repmat(1:N,Kr,1);              % [Kr x N] index into ROIs
        end
        

        % ABC 2017-05-18: BUG fix: should be N (usually Kr>N for tagging)
        [~, ~, lp_out] = GMM_MixHierEM(trainMix, tmpK, virtualSamples, iterations);
        %[reduced_out] = GMM_MixHierEM(trainMix, Kr, virtualSamples, iterations);
        % post_out is the position,
        %lp_out is like the hat_z_ij


        % map from HMM emission to group ROI
        emit_map_ll = cell(1,Kb);
        for i=1:numbc
          emit_map_ll{gmms_hmm_index(i)}(:,gmms_hmm_emit_index(i)) = lp_out(:,i);
        end

        % assign base HMMS to reduced HMM using emissions
        hmm_map = zeros(1,Kb);
        hmm_map_emit = cell(1,Kb);
        hat_z = zeros(Kb,Kr);

        for i=1:Kb
          % compute sum of LL between group ROIs and HMM emissions
          emit_LL = zeros(1,Kr);
          for j=1:Kr
            emit_LL(j) = sum(sum(emit_map_ll{i}(useGauss(j,:),:)));
            % the emission of jth HMM in reduced model consist of useGauss(j,:)
          end
          % pick largest sum LL as the group
          [~, myj] = max(emit_LL);
          % map the base ROIs to group ROIs.
          [~, myii] = max(emit_map_ll{i}(useGauss(myj,:),:),[], 1);   %?
          hmm_map(i) = myj;
          hat_z(i,myj) = 1;
          hmm_map_emit{i} = myii;
        end


        % build new transition matrices
        newA = cell(1,Kr);
        newprior = cell(1,Kr);
        newN = zeros(1,Kr);   
        for j=1:Kr
          newA{j} = zeros(N,N);
          newprior{j} = zeros(N,1);
          newN(j) = 0;
        end
        for i=1:Kb      
          mymap = hmm_map_emit{i}; % get the ROI map from base to reduced
          myj = hmm_map(i);      % get the corresponding group
          % for each base hidden state
          for k=1:length(mymap)
            nk = mymap(k);
            % add probability to corresponding reduced hidden state
            newprior{myj}(nk) = newprior{myj}(nk) + h3m_b.hmm{i}.prior(k);
            for l=1:length(mymap)
              nl = mymap(l);
              newA{myj}(nk,nl) = newA{myj}(nk,nl) + h3m_b.hmm{i}.A(k,l);
            end
          end
          newN(myj) = newN(myj)+1;
        end
        % normalize
        for j=1:Kr
          newprior{j} = newprior{j} / sum(newprior{j});
          newA{j} = bsxfun(@rdivide, newA{j}, sum(newA{j}+ 1e-50,2));
        end
        
        % start random
    
        wb      = h3m_b.omega; %size =1*Kb
        tilde_N_k =  virtualSamples *wb'; %Kb*1 

        % seperate the z by let the max possible z=1/2, and other average the other 1/2 
        [maxz,idz] = max(hat_z,[],2);
        hat_zz = zeros(size(hat_z));
        for i = 1:Kb
            hat_zz(i,:)=maxz(i)/(2*(Kr-1));
            hat_zz(i,idz(i))= maxz(i)/2;
        end
        z_N = hat_zz.*repmat(tilde_N_k,1,Kr);

        % alpha
        N_ll = sum(z_N,1);        
        alpha1 = N_ll + alpha0;
        h3m.alpha = alpha1;

        Nu_1           = cell(Kr,Kb);
        Sum_t_nu       = cell(Kr,Kb);
        Sum_xi         = cell(Kr,Kb);

        % find the max possible cluster
        idmax_z = sum(hat_z,1)==max(sum(hat_z,1));
        for i = 1: Kr    
            if ~any(newA{1,i}(:))
                newA{1,i} = newA{1,idmax_z};
                newprior{1,i} = newprior{1,idmax_z};
            end
        end

        for i = 1: Kr
            
            for j= 1:Kb

                Sb = size(h3m_b.hmm{j}.emit,2);
                
                Theta_1 =  newprior{i} * ones(1,Sb);

                Theta =zeros(N,N,Sb,T);

                for t=2:T
                    for beta = 1: Sb      
                        Theta(:,:,beta,t) = newA{i};
                    end
                end

                %get pars of jth hmm
                pib = h3m_b.hmm{j}.prior;   %  column vect
                Ab = h3m_b.hmm{j}.A;

                nu = (ones(N,1) * pib') .* Theta_1;
                sum_nu_1 = sum(nu,2)';    %hat v_1

                sum_t_nu = nu;
                % CACHE: sum_t sum_gamma xi(rho,sigma,gamma,t)
                sum_t_sum_g_xi = zeros(N,N); % N2 by N2 (indexed by rho and sigma) hat_xi

                for t = 2 : T
                    % compute the inner part of the update of xi (does not depend on sigma)
                    foo = nu * Ab; % indexed by rho gamma

                    xi_foo_all = bsxfun(@times, reshape(foo,[N 1 Sb]), Theta(:,:,:,t));

                    % sum_t_sum_g_xi(:,sigma) = sum_t_sum_g_xi(:,sigma) + sum(xi_foo,2);
                    sum_t_sum_g_xi = sum_t_sum_g_xi + sum(xi_foo_all, 3);

                    % nu(sigma,:) = ones(1,N2) * xi_foo;
                    nu = reshape(sum(xi_foo_all,1), [N Sb]);
  
                     % CACHE: in the sum_t nu_t(sigma, gamma)
                    sum_t_nu = sum_t_nu + nu;    
                end

                sum_xi = sum_t_sum_g_xi;
                Nu_1{i,j} = sum_nu_1;
                Sum_xi{i,j} = sum_xi;
                Sum_t_nu{i,j} = sum_t_nu;
           end
        end

        up_eta = zeros(Sr,Kr);
        up_eps = zeros(Sr,Sr,Kr);
        Nl_j     = zeros(Sr,Kr);       
        y_bar    = zeros(Sr,dim,Kr);
        update_S = zeros(dim,dim,Sr,Kr);
        switch vbhemopt.emit.covar_type
            case 'diag'
                update_C = zeros(Sr,dim,Kr);        
        end

        for j = 1 : Kr

               % loop all the components of the base mixture
            for i = 1 : Kb

                Sb = size(h3m_b.hmm{i}.emit,2);        
                mubl = zeros(dim,Sb);
                Mubl = zeros(Sb,dim);
                   % for each emission density, convert to H3M format
                for h = 1:Sb            
                    mubl(:,h) =h3m_b.hmm{i}.emit{h}.centres';
                    Mubl(h,:) = h3m_b.hmm{i}.emit{h}.covars;
                end  

                nu     = Nu_1{j,i};              % this is a 1 by N2 vector
                xi     = Sum_xi{j,i};            % this is a N2 by N2 matrix (from - to)
                hat_nu = Sum_t_nu{j,i};

                up_eta(:,j) = up_eta(:,j) + z_N(i,j) * nu';   % up eta
                up_eps(:,:,j) = up_eps(:,:,j) + z_N(i,j) * xi;     % up epsilon
                Nl_j(:,j)  = Nl_j(:,j) + z_N(i,j) * sum(hat_nu,2)+ 1e-50;

                y_bar(:,:,j) = y_bar(:,:,j) + z_N(i,j) * (hat_nu*mubl');  %row

                switch vbhemopt.emit.covar_type
                    case 'diag'  
                        update_C(:,:,j)  =update_C(:,:,j) + z_N(i,j) * (hat_nu*Mubl);
                end   
                for k = 1:Sr
                    for h = 1: Sb
                        update_S(:,:,k,j) = update_S(:,:,k,j) + z_N(i,j) * (hat_nu(k,h)*(mubl(:,h)*mubl(:,h)'));
                    end
                end    
            end 

            y_bar(:,:,j) = y_bar(:,:,j)./(Nl_j(:,j)*ones(1,dim));
            for i=1:Sr
                update_S(:,:,i,j) =  update_S(:,:,i,j)/Nl_j(i,j) - y_bar(i,:,j)'*y_bar(i,:,j);
            end
            switch vbhemopt.emit.covar_type  
                case 'diag'
                    update_C(:,:,j)  =update_C(:,:,j) ./(Nl_j(:,j)*ones(1,dim));
            end                  
        end
        N_Eta0 = up_eta;
        N_Eps0 = up_eps;

        eta = eta0 + N_Eta0;    % Sr * Kr
        lambda = lambda0 + Nl_j;   % Sr * Kr
        v = v0 + Nl_j + 1;

        for i = 1:Kr
          h3m.hmms{i}.eta = eta(:,i);
          epsilon = N_Eps0(:,:,i) + epsilon0;
          h3m.hmms{i}.epsilon = epsilon;    

          h3m.hmms{i}.lambda = lambda(:,i);
          h3m.hmms{i}.v = v(:,i);

          m = zeros(dim,Sr);
          W = zeros(dim,dim,Sr);
          Cov = zeros(dim,dim,Sr);
          h3m.hmms{i}.emit = {};
          for k = 1:Sr
            m(:,k) = (lambda0.*m0 + Nl_j(k,i).*y_bar(k,:,i)')./lambda(k,i);
            h3m.hmms{i}.emit{k}.centres = m(:,k)';
            mult1 = lambda0.*Nl_j(k,i)/(lambda0 +Nl_j(k,i));
            diff3 = y_bar(k,:,i) - m0';
            W(:,:,k) = inv(W0inv + Nl_j(k,i)*(diag(update_C(k,:,i))  + ...
                update_S(:,:,k,i)  ) + mult1*diff3'*diff3  );
            h3m.hmms{i}.emit{k}.covars = diag(W(:,:,k));

            h3m.hmms{i}.emit{k}.v = v(k,i);
            h3m.hmms{i}.emit{k}.lambda = lambda(k,i);

          end

          % calculate covariance matrices
          for k = 1:Sr
              Cov(:,:,k) = inv(W(:,:,k))/(v(k)-dim-1);
          end
          h3m.hmms{i}.m = m; 
          h3m.hmms{i}.W = W;
          h3m.hmms{i}.Cov = Cov;
        end

        h3m.hmm = h3m.hmms;
        h3m_r = h3m;


        
    case 'con-wtkmeans'      %' 
        
        
        % make the initialization more consistent
        % at first, provide a initialization which is produced by wtkeamns
        % with largst K and S, then, cluster this centers to small K and S
        
        
        cluster_center = vbhemopt.cluster_center;
        cluster_idx = vbhemopt.cluster_idx;
        cluster_centers = vbhemopt.cluster_centers;
        
        weight_idx = hist(cluster_idx,unique(cluster_idx))/length(cluster_idx);
        % give a start seed for the kmeans
        it_max =100;
        rng(vbhemopt.wtseed, 'twister');
        [~,init_center] = kmeans(cluster_center',Kr,'Display','off','Replicates',1);    %'input is n*dim, output is k*dim
        % use Kmeans results as initial for weighted kmeans
        [cluster_idx2,~] = my_weighted_kmeans(Kr,it_max,cluster_center,weight_idx,init_center');
        cluster_centers_new = cell(Kr,1);
        empty_ind = zeros(Kr,1);

        if length(unique(cluster_idx2)) < (Kr/2)
            
            cluster_centers_new = cellfun(@(x) x(:,1:Sr),cluster_centers(1:Kr), 'uniformoutput', false);

        else
            
        for i = 1:Kr
            if (sum(cluster_idx2==i)==0)
                         
                cluster_centers_new{i} = [];
                empty_ind(i) = 1;
                continue
                    
            else
            temp = cluster_centers(cluster_idx2==i);
            
            temp2 = cat(2,temp{:});  %dim*n
            if size(temp2,2)<Sr
                init_center2 = temp2';  %n*dim
            else
                [~,init_center2] = kmeans(temp2',Sr,'Display','off','Replicates',1);    %'input is n*dim, output is k*dim
            end
            if size(init_center2,1)<Sr
                init_center2(end+1:Sr,:) = init_center2(1,:);
            end
            cluster_centers_new{i,1} = init_center2';

            end                
        end        
        end
        

        cluster_centers = cluster_centers_new;
        
        first_unempty = find(empty_ind~=1, 1, 'first');
        if sum(empty_ind)~=0
          for i = 1:Kr       
              if empty_ind(i)
                  cluster_centers{i} = cluster_centers{first_unempty};
              end
          end
        end
        

        virtualSamples = Nv; 
        Nv = ones(1,Kr)/Kr* virtualSamples;
        h3m = h3m;
        
        h3m.hmms = {};
    
        for j=1:Kr
            
          [prior, A] = makeAprior(N, vbhemopt.initopt.mode);

          h3m.hmms{j}.eta   = prior * Nv(j) + eta0 ;  %initial uniform P(pi_rho1) = 1/Kr 
          h3m.hmms{j}.epsilon  = A *Nv(j) + epsilon0; 
          Nsj_rho = Nv(j)*ones(1,Sr)/Sr;                    
          h3m.hmms{j}.lambda = lambda0 + Nsj_rho';
          h3m.hmms{j}.v = v0 + Nsj_rho' + 1;
          h3m.hmms{j}.emit = {};
          
          cluster_center = cluster_centers{j};
          
          if size(cluster_center,2)<Sr
              keyboard
          end
          for n=1:Sr

              h3m.hmms{j}.m(:,n)= cluster_center(:,n);
              h3m.hmms{j}.emit{n}.centres = cluster_center(:,n)';
              h3m.hmms{j}.emit{n}.v = h3m.hmms{j}.v(n);
              h3m.hmms{j}.emit{n}.lambda = h3m.hmms{j}.lambda(n);

              switch(h3m_b.hmm{1}.emit{1}.covar_type)
              case 'diag'
                 h3m.hmms{j}.W(:,:,n)= inv((h3m.hmms{j}.v(n)-dim-1) *diag((h3m_b.hmm{1, 1}.emit{1, 1}.covars)));
                 h3m.hmms{j}.emit{n}.covars = diag(h3m.hmms{j}.W(:,:,n));
                 h3m.hmms{j}.Cov(:,:,n)= h3m_b.hmm{1, 1}.emit{1, 1}.covars;
              case 'full'
                  h3m.hmms{j}.W(:,:,n)= inv(W0inv  + Nsj_rho(n) *((h3m.hmms{j}.emit{n}.covars)));
              end    
          end      
        end
    
        h3m.alpha = alpha0 + Nv;
        h3m.hmm = h3m.hmms;
        h3m_r =h3m;
        
        
 
    case 'random'
        
        % given a random assignment, using the means of the hmm that assigned to
        % the j hmm to fit the emmition gmm. trans and ini is average

        h3m_b_hmm_emit_cache = make_h3m_b_cache(h3m_b);
        %random z
        post = randi(Kr,[1,Kb]);
        % get an assignment that assign data to each component

        % get rid of all data assign to the same component
        while 1
            if Kr~=1
                if Kr ~= length(unique(post))  
                    post = randi(Kr,[1,Kb]);
                else
                    break
                end
            else
                break
            end
        end


        Rad_z = zeros(Kb,Kr);
        mix_all = cell(1,Kr);
        for j = 1:Kr
            mix = struct;
            mu_temp = [];
            ind_b= find(post==j);
            for i = 1:length(ind_b)
                mu_temp = [mu_temp,reshape(h3m_b_hmm_emit_cache{i}.gmmBcentres,[dim,length(h3m_b_hmm_emit_cache{i}.gmmBpriors)])];
            end

            mu_temp = mu_temp';   %n*dim'
            Nd = size(mu_temp,1);

            if (Sr==1)
                % check if just one component, which is a Gaussian
                mix.PComponents = [1];
                mix.mu    = mean(mu_temp);
                mix.Sigma = cov(mu_temp);

            elseif (Nd<=Sr)
              % check if enough data to learn a GMM
              % if not, then
              %   select point as the means and iid variance
              Nextra = Sr-Nd;
              tmp = [ones(1,Nd), 0.000001*ones(1,Nextra)];      
              mix.PComponents = tmp / sum(tmp);      
              mix.mu    = [mu_temp; zeros(Nextra, dim)];
              tmp = var(mu_temp);
              mix.Sigma = mean(tmp)*repmat(eye(dim), [1 1 Sr]);

            else
              % enough data, so run GMM
              try
                warning('off', 'stats:gmdistribution:FailedToConverge');

                if ~isempty(vbhemopt.random_gmm_opt)
                  if (VERBOSE_MODE >= 3)
                    fprintf('random gmm init: passing random_gmm_opt = ');
                    vbhemopt.random_gmm_opt
                  end
                end

                % use the Matlab fit function to find the GMM components
                % 2017-01-21 - ABC - added ability to pass other options (for Antoine)
                % 2018-11-26 - v0.74 - to be consistent with previous versions, specify 'Start' as 'randSample' (which was the default previously)
                mix = gmdistribution.fit(mu_temp,Sr, vbhemopt.random_gmm_opt{:}, ...
                  'Options', struct('TolFun', 1e-5), 'Start', 'randSample');

              catch ME
                if strcmp(ME.identifier, 'stats:gmdistribution:IllCondCovIter')
                  if (VERBOSE_MODE >= 2)
                    fprintf('using shared covariance');
                  end
                  % 2018-02-01: added regularize option for stability.
                  mix = gmdistribution.fit(mu_temp, Sr, 'SharedCov', true, 'Regularize', 1e-10, ...
                    'Options', struct('TolFun', 1e-5), 'Start', 'randSample');
                  %mix = gmdistribution.fit(data,K, 'SharedCov', true, 'Options', struct('TolFun', 1e-5));

                  % SharedCov -> SharedCovariance

                else

                  % otherwise use our built-in function (if Stats toolbox is not available)
                  if (VERBOSE_MODE >= 2)
                    fprintf('using built-in GMM');                        
                    warning('vbhmm_init:nostatstoolbox', 'Stats toolbox is not available -- using our own GMM code. Results may be different than when using the Stats toolbox.');
                    warning('off', 'vbhmm_init:nostatstoolbox');
                  end
                  gmmopt.cvmode = 'full';
                  gmmopt.initmode = 'random';
                  gmmopt.verbose = 0;
                  gmmmix = gmm_learn(mu_temp', Sr, gmmopt);
                  mix.PComponents = gmmmix.pi(:)';
                  mix.mu = cat(2, gmmmix.mu{:})';
                  mix.Sigma = cat(3, gmmmix.cv{:});
                end
              end
            end

            mix_all{j} = mix;
            Rad_z(ind_b,j) = 1; 
        end

        N_i = Nv *h3m_b.omega;
        N_ij = Rad_z.*N_i';
        N_j = sum(N_ij);
        N_j2 = N_j.*repmat(1/Sr,Sr,1);

        h3m.alpha = alpha0 + N_j;  % row

        for j = 1:Kr

            mixj = mix_all{j};        
            hmm.eta   =  eta0 + N_j2(:,j); %column

            for k = 1:Sr           
                hmm.epsilon(k,:) = epsilon0 + N_j2(:,j)';
            end

            % Update gmm
            Nj_rho = N_j(j)*mixj.PComponents';
            hmm.lambda = lambda0 + Nj_rho;
            hmm.v = v0 + Nj_rho +1;

            ybar = mixj.mu';
            hmm.m = ((lambda0*m0)*ones(1,Sr) + (ones(dim,1)*Nj_rho').*ybar)./(ones(dim,1)*hmm.lambda');

            if size(mixj.Sigma,3) == Sr
                Sig = mixj.Sigma;
            elseif (size(mixj.Sigma,3) == 1)
                % handle shared covariance
                Sig = repmat(mixj.Sigma, [1 1 Sr]);
            end


           % for diagonal covarainces, one of the dimensions will have size 1
           if ((size(Sig,1) == 1) || (size(Sig,2)==1)) && (dim > 1)

               oldS = Sig;
               Sig = zeros(dim, dim, K);
               for k=1:Sr               
                   Sig(:,:,k) = diag(oldS(:,:,k)); % make the full covariance
               end
           end  


           hmm.W = zeros(dim,dim,Sr);

           for k = 1:Sr           
               mult1 = lambda0.*Nj_rho(k)/(lambda0 + Nj_rho(k));
               diff3 = ybar(:,k) - m0;
               hmm.W(:,:,k) = inv(W0inv + Nj_rho(k)*Sig(:,:,k) + mult1*diff3*diff3');

           end

             h3m.hmm{j} = hmm;

        end

        h3m.LogL = -inf;
        h3m_r = h3m;
    
    
          %%  case
       case 'vbh3m'   %given an H3M , using for compare with vhem, using the same initial h3m.
        
        vbh3m = vbhemopt.given_h3m;%vbhemopt.vbh3m{Kr, S};
        if ~isfield(vbh3m,'hmm')
            vbh3m.hmm = vbh3m.hmms;
        end
         
        Weight = vbh3m.group_size./sum(vbh3m.group_size);
        for j=1:Kr    
            NLrj = (Nv*Weight(j));
            S = length(vbh3m.hmm{j}.prior);
            NLrjn =[];
          for n = 1:S    % state
            NLrjn(n) = NLrj*vbh3m.hmm{j}.prior(n);
           % h3m.hmm{j}.emit{n}.centres = vbh3m.hmm{j}.pdf{n}.mean;
            h3m.hmm{j}.emit{n}.m = vbh3m.hmm{j}.pdf{n}.mean;
            h3m.hmm{j}.emit{n}.W = inv( (v0 + NLrjn(n) -dim)*vbh3m.hmm{j}.pdf{n}.cov); 
         %   h3m.hmm{j}.emit{n}.covars = inv( (v0 + NLrjn(n) -dim)*vbh3m.hmm{j}.pdf{n}.cov);  %W not cov, here cov is a matrix
%             h3m.hmm{j}.m(:,n)= vbh3m.hmm{j}.pdf{n}.mean';
%             h3m.hmm{j}.W(:,:,n)= inv( (v0 + NLrjn(n) -dim)*vbh3m.hmm{j}.pdf{n}.cov); 
            h3m.hmm{j}.emit{n}.lambda = lambda0 + NLrjn(n) ;
            h3m.hmm{j}.emit{n}.v = v0 + NLrjn(n) +1;
          end
          %[h3m.hmms{j}.prior, h3m.hmms{j}.A] = makeAprior(N, vbhemopt.initopt.mode);
          h3m.hmm{j}.eta   = NLrjn' + eta0 ;  %initial uniform P(pi_rho1) = 1/Kr 
          h3m.hmm{j}.epsilon = (vbh3m.hmm{j}.trans .* repmat(NLrjn',1,S)) + epsilon0;   
%           h3m.hmm{j}.lambda = lambda0 + NLrjn';
%           h3m.hmm{j}.v = v0 + NLrjn'+1;

        end
        
%          omega = rand(Kr,1);
% %   omega = ones(Kr,1)/Kr;
             omega =Weight';
%         omega = omega/sum(omega);
%         %omega = ones(1,Kr)/Kr;
        h3m.LogL = -inf;
        h3m.alpha = alpha0 + omega'*Nv;
        h3m_r =h3m;
        
    case 'given_h3m'
    
        vbh3m = vbhemopt.given_h3m;
        h3m = h3m;
        h3m.hmms = vbh3m.hmm(1:Kr);
    
        for j=1:Kr

          h3m.hmms{j}.eta   = vbh3m.hmm{1, j}.eta(1:Sr);  %initial uniform P(pi_rho1) = 1/Kr 
          h3m.hmms{j}.epsilon  = vbh3m.hmm{1, j}.epsilon(1:Sr,1:Sr);                    
          h3m.hmms{j}.lambda = vbh3m.hmm{1, j}.lambda(1:Sr);
          h3m.hmms{j}.v = vbh3m.hmm{1, j}.v(1:Sr);
          h3m.hmms{j}.emit = vbh3m.hmm{1, j}.emit(1:Sr);
          h3m.hmms{j}.m = vbh3m.hmm{1, j}.m(:,1:Sr);
          h3m.hmms{j}.W = vbh3m.hmm{1, j}.W(:,:,1:Sr);
          h3m.hmms{j}.Cov = vbh3m.hmm{1, j}.W(:,:,1:Sr);
        end
        
        
        h3m.alpha = vbh3m.alpha(1:Kr);
        h3m.hmm = h3m.hmms;
        h3m_r =h3m;
        
        
end


function [p_out] = makeGMMweights(p, A, mode)
p = p(:)';
  
if (any(mode=='0'))
  for t = 1 : 50
    p = p * A;
  end
  p_out = p;
  
elseif (any(mode=='3'))
  T = 3;
  p_out = p;
  for t=2:T
    p = p * A;
    p_out = p_out + p;
  end  
  p_out = p_out/T;

elseif (any(mode=='x'))
  p_out = ones(size(p))/size(p,2);

else
  error('unknown mode');
end




%% make the prior and A
function [prior, A] = makeAprior(N, mode)
if any(mode == 'u')
  prior = ones(N,1)/N;
  A     = ones(N,N)/N;
elseif any(mode == 'r')
  prior = rand(N,1);
  prior = prior/sum(prior);
  A = rand(N);
  A = A ./ repmat(sum(A,2),1,N);
else
  error('unknown mode');
end

