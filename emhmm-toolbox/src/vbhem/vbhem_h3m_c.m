function [h3m_new] = vbhem_h3m_c(h3m_b, vbhemopt) 
% this function to reduce the h3m_b to a h3m_new with K components and S
% states
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-07-29
% Antoni B. Chan, Janet H. Hsiao, Lan Hui
% City University of Hong Kong, University of Hong Kong

% 2021-07-29: v0.80 - initial version


VERBOSE_MODE = vbhemopt.verbose;
Verbose_Prefix =vbhemopt.verbose_prefix;
% check for estimating hyps automatically
if iscell(vbhemopt.learn_hyps) || (vbhemopt.learn_hyps == 1)
  do_learn_hyps = 1;
else
  do_learn_hyps = 0;
end

% 
K = vbhemopt.K;
S = vbhemopt.S;

numits   =  vbhemopt.trials;

% if K==1
%     numits=1;
% end
h3m_news = cell(1,numits);
LLall    = zeros(1,numits);

%% do multiple trials and keep the best
fprintf('VBHEM Trial: ');
parfor it = 1 : numits            %% run the t times , not really time
  vbhemopt2 = vbhemopt;
  % set the seed here
  if ~isempty(vbhemopt2.seed)
    if (VERBOSE_MODE > 2)
      fprintf('%s+ setting random seed to %d\n', vbhemopt.verbose_prefix, vbhemopt.seed);
    elseif (VERBOSE_MODE == 1)
      % fprintf('(seed=%d)', vbhemopt2.seed+it);
      rng(vbhemopt2.seed+it, 'twister');
      
    end
  end
  
  h3m_r = [];
  
  % initialization
  switch vbhemopt2.initmode
    
    %     case 'given_h3m'
    %        %vbhemopt2 = vbhemopt;
    %        vbhemopt2.given_h3m = vbhemopt.given_h3m{it};
    %        h3m_r = vbhemhmm_init(h3m_b,vbhemopt2);
    
    
    case {'random','baseem', 'gmmNew', 'gmmNew2','gmmNew2_r','wtkmeans','vbh3m'}
      
      if   strcmp(vbhemopt2.initmode, 'wtkmeans')
        vbhemopt2.wtseed = vbhemopt.seed + it;
      end
      if   isfield(vbhemopt2, 'all_rands')
        vbhemopt2.prerand = vbhemopt.all_rands{it};
      end
      
      [h3m_r,vbhemopt2] = vbhemhmm_init(h3m_b,vbhemopt2);
      
    otherwise
      error('bad init')
  end
  
  
  h3m_news{it} = vbhem_h3m_c_step_fc(h3m_r,h3m_b,vbhemopt2);
  
  LLall(it) = h3m_news{it}.LL;
  
  if (VERBOSE_MODE == 1)
    fprintf('.');
    if mod(it,10)==0
      fprintf('%d', it);
    end
  elseif (VERBOSE_MODE > 1)
    fprintf('Trial %d: LL=%g\n', it, LLall(it));
  end
  
end

show_LL_plot = 0;
% get unique solutions

if (vbhemopt.keep_suboptimal_hmms || do_learn_hyps)
% if percentage change is larger than twice convergence criteria.
% 2*   because two models could be on either side of the true model.
% 10*  since convergence might not happen completely
    diffthresh = 2*vbhemopt.minDiff*10;

    [unique_h3m_it, VLL, legs] = uniqueLL(LLall, diffthresh, show_LL_plot);

    unique_h3m    = h3m_news(unique_h3m_it);
    unique_h3m_LL = LLall(unique_h3m_it);
end

%% learn hyps after random initializations


if (do_learn_hyps)
    
    show_hyp_plot = 0;

    % save old values
    LLall_old    = LLall;
    h3m_news_old = h3m_news;

    % invalidate
    LLall = nan*LLall_old;
    h3m_news = cell(size(h3m_news_old ));
    h3m_news1 = cell(length(unique_h3m),1);
    
    parfor q=1:length(unique_h3m)
      
      my_it  = unique_h3m_it(q);
      my_h3m = unique_h3m{q};
      my_LL  = unique_h3m_LL(q);

      %fprintf('\noptimizing trial %d: %g', my_it, my_LL);
      if (VERBOSE_MODE >= 1)
        fprintf('%s(K=%d) optimizing trial %d: %g', Verbose_Prefix, K, my_it, my_LL);
      end

      % estimate hyperparameters
      vbhemopt2 = vbhemopt;
      vbhemopt2.initmode = 'inith3m';
      vbhemopt2.inithmm  = my_h3m;
%       h3m_news{my_it} = vbhem_h3m_c_hyp(h3m_b, vbhemopt2);
%       LLall(my_it)  = h3m_news{my_it}.LL;

      h3m_news1{q} = vbhem_h3m_c_hyp(h3m_b, vbhemopt2);
      LLall_1(q)  = h3m_news1{q}.LL;


      % visualization
      if (show_LL_plot)
        figure(VLL)
        legs4 = plot(my_it, LLall(my_it), 'rx');
        legend([legs, legs4], {'original LL', 'selected trial', 'optimized LL', 'ignore bounds'});
        drawnow
      end
    end

    if (show_hyp_plot)
     % show models
     faceimg ='face.jpg';
     figure         
     NN = length(unique_h3m);
     for i=1:NN
       myit = unique_h3m_it(i);
       subplot(2,NN,i)
       vbhmm_plot_compact(unique_h3m{i}, faceimg);
       title(sprintf('trial %d\norig LL=%g', myit, unique_h3m{i}.LL));
       subplot(2,NN,i+NN)
       vbhmm_plot_compact(h3m_news{myit}, faceimg);
       title(sprintf('opt LL=%g', h3m_news{myit}.LL));
     end
     drawnow
    end
    
    h3m_news(unique_h3m_it) = h3m_news1(:,1);
    LLall(unique_h3m_it) = LLall_1;

end % end do_learn_hyps


%% choose the best
[maxLL,maxind] = max(LLall);
h3m_new = h3m_news{maxind};
t_best = maxind;
%it_best = h3m_new.iter;
LL_best = h3m_new.LL;  
if (VERBOSE_MODE >= 1)  
fprintf('Best run is %d: LL=%g\n',t_best,LL_best)
end

if (do_learn_hyps)
    
   % check degenerate models and save (DEBUGGING)
   if any(abs(LLall_old./LLall)>10)
       if (VERBOSE_MODE >= 3)
         foo=tempname('.');
         warning('degenerate models detected: saving to %s', foo);
         save([foo '.mat'])
       else
         warning('degenerate models detected');
       end
   end

    % output final params
    if (VERBOSE_MODE >= 1)
      fprintf('  ');          
      vbhmm_print_hyps({h3m_new.learn_hyps.hypinfo.optname}, h3m_new.learn_hyps.vbopt);
      fprintf('\n');
    end

%     % choose the best among the random initialization
%     if (vbhemopt.keep_best_random_trial)
%         [maxLL_old,maxind_old] = max(LLall_old);
%         tmph3m = h3m_news_old{maxind_old};
%         if ~isempty(vbhemopt.sortclusters)
%             for j = 1:tmph3m.K 
%                 tmphmm = tmph3m.hmm{j};
%                 tmphmm = vbhmm_standardize(tmphmm, vbhemopt.sortclusters);
%                 tmph3m.hmm{j} = tmphmm;
%             end
%         end
%         h3m_new.learn_hyps.hmm_best_random_trial = tmph3m;
%     end
end

if (vbhemopt.keep_suboptimal_hmms)
h3m_new.suboptimal_hmms = unique_h3m;
end





