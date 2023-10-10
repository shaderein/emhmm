function [h3m_new] = vbhem_coh3m_c(h3m_bs, vbhemopt) 
% this function to reduce the h3m_b to a h3m_new with K components and S
% states (internal function)
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-07-29
% Antoni B. Chan, Janet H. Hsiao, Lan Hui
% City University of Hong Kong, University of Hong Kong

% 2021-07-29: v0.80 - Hui - initial version

VERBOSE_MODE = vbhemopt.verbose;
Verbose_Prefix =vbhemopt.verbose_prefix;

if (vbhemopt.K == vbhemopt.subjectsize) || (vbhemopt.K == 0)
    h3m_new = h3m_bs;
    return
end

% check for estimating hyps automatically
if iscell(vbhemopt.learn_hyps) || (vbhemopt.learn_hyps == 1)
  do_learn_hyps = 1;
else
  do_learn_hyps = 0;
end

K = vbhemopt.K;

Kb = vbhemopt.subjectsize;
Nstimuli = vbhemopt.stimulisize;
% ABC 2020-08-06: BUG FIX v0.xx: remove empty HMMs for performing initialization
% this removes the association between subjects between stimuli, but is
% fine for initialization.
h3m_bs_no_empty = h3m_bs;
validinds = cell(Nstimuli,1);
fprintf('init check: ignore missing hmms: ');
for i=1:Nstimuli 
  myinds = find(h3m_bs(i).omega > 0);
  validinds{i} = myinds;
  fprintf(' %d', Kb-length(myinds));
  
  if (length(myinds) < Kb)
    h3m_bs_no_empty(i).K = length(myinds);
    h3m_bs_no_empty(i).omega = h3m_bs_no_empty(i).omega(myinds);
    h3m_bs_no_empty(i).omega = h3m_bs_no_empty(i).omega / sum(h3m_bs_no_empty(i).omega);
    h3m_bs_no_empty(i).hmm = h3m_bs_no_empty(i).hmm(myinds);
  end
end
fprintf('\n');


numits   =  vbhemopt.trials;
h3m_news = cell(1,numits);
LLall    = nan*zeros(1,numits);
%All_inihmms = cell(numits,1);

% do multiple trials of co-clustering
fprintf('VBHEM Coclustering Trial: ');
parfor it = 1 : numits            %% run the t times , not really time
  
  if ~isempty(vbhemopt.seed)
    %fprintf('+ seed %d\n', mopt.seed+t);
    seed_thisit = vbhemopt.seed+it;
    rng(seed_thisit, 'twister');
  end
  
  
  %  2020-08-06: ignore empty HMMs
  h3m_bs_pointer = h3m_bs_no_empty;
  
  h3ms = [];
  
  % initialization
  switch vbhemopt.initmode
    
    %     case 'given_h3m'
    %        vbhemopt2 = vbhemopt;
    %        vbhemopt2.given_h3m = vbhemopt.given_h3m{it};
    %        h3m_r = vbhemhmm_init(h3m_bs,vbhemopt2);
    %
    %     case {'random','baseem', 'gmmNew', 'gmmNew2','gmmNew2_r','wtkmeans'}
    %
    %         try
    %             % use less data (400) for initialization
    %             Kconstant = 400;
    %             h3ms = [];
    %             moptpointer = vbhemopt;
    %             moptpointer.keep_suboptimal_hmms = 0;
    %
    %             if Kb > Kconstant
    %                 inds = randi(Kb-Kconstant);
    %                 new_idx = inds:(inds+Kconstant-1);
    %                 %initialize all stimuli
    %                 for i=1:vbhemopt.stimulisize
    %                     h3m_bs_pointer(i).K = Kconstant;
    %                     h3m_bs_pointer(i).omega = h3m_bs_pointer(i).omega(new_idx);
    %                     h3m_bs_pointer(i).omega = h3m_bs_pointer(i).omega / sum(h3m_bs_pointer(i).omega);
    %                     h3m_bs_pointer(i).hmm = h3m_bs_pointer(i).hmm(new_idx);
    %                     if length(vbhemopt.S)>1
    %                         moptpointer.S = vbhemopt.S(i);
    %                     end
    %                     h3m = vbhemhmm_init(h3m_bs_pointer(i),moptpointer);
    %                     h3ms = [h3ms; h3m];
    %                 end
    %
    %
    %             else
    %
    %                 %initialize all stimuli
    %                 for i = 1:vbhemopt.stimulisize
    %                     if length(vbhemopt.S)>1
    %                         moptpointer.S = vbhemopt.S(i);
    %                     end
    %                %     moptpointer.given_h3m = vbhemopt.given_h3m{K}{i};
    %                     h3m = vbhemhmm_init(h3m_bs_pointer(i),moptpointer);
    %                     h3ms = [h3ms; h3m];
    %                 end
    %
    %             end
    %
    %             if (isfield(vbhemopt,'keep_ini_hmms')&&(vbhemopt.keep_ini_hmms))
    %               All_inihmms{it} = h3ms;
    %             end
    %
    %         catch me
    %             % there is some bug in arkmeans -- need to change it
    %             % anyways, invalidate this trial
    %             h3ms = -1;
    %           %  keyboard
    %         end
    %
    case 'vbh3m'
      
      %h3ms = struct;
      h3ms =[];
      moptpointer = vbhemopt;
      moptpointer.keep_suboptimal_hmms = 0;
      ini_cogroup_hmms = vbhemopt.VHintial{it};
      
      %initialize all stimuli
      for i = 1:Nstimuli
        if length(vbhemopt.S)>1
          moptpointer.S = vbhemopt.S(i);
        end
        moptpointer.given_h3m = ini_cogroup_hmms{i};
        h3m = vbhemhmm_init(h3m_bs_pointer(i),moptpointer);
        h3ms = [h3ms; h3m];
      end
      
      
  end
  
  % initialize h3m_new (for parfor loop)
  h3m_new = struct;
  h3m_new.LL = -inf;
  
  if isnumeric(h3ms) && (h3ms == -1)
    % bad initialization, so skip this HEM run
    h3m_new = struct;
    h3m_new.LL = -inf;
    
  else
    
    %% run HEM using h3m as the initialization
    % the new reduced h3m is two-dimentional and different from the
    % orginal H3M format
    try
      h3m_new = vbhem_coh3m_c_step_fc(h3ms, h3m_bs, vbhemopt);
            
    catch me
      h3m_new = struct;
      h3m_new.LL = -inf;
    end
    
  end
  
  % save this iteration
  h3m_news{it} = h3m_new;
  LLall(it) = h3m_new.LL;
  
  if (VERBOSE_MODE == 1)
    fprintf('.');
    if mod(it,10)==0
      fprintf('%d', it);
    end
  elseif (VERBOSE_MODE > 1)
    fprintf('Trial %d: LL=%g\n', it, LLall(it));
  end
  
end

clear h3m_bs_pointer
clear h3m_bs_no_empty
clear moptpointer

show_LL_plot = 0;

%% get unique solutions
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
  
  %  show_hyp_plot = 0;
  
  % save old values
  LLall_old    = LLall;
  % h3m_news_old = h3m_news;
  
  
  % invalidate
  sizeh3m_news = size(h3m_news);
  clear h3m_news
  LLall = nan*LLall;
  h3m_news1 = cell(length(unique_h3m),1);
  
  
  % for each iteration w/ unique h3m
  % parfor
  parfor q = 1:length(unique_h3m)
    
    my_it  = unique_h3m_it(q);
    my_h3m = unique_h3m{q};
    my_LL  = unique_h3m_LL(q);
    
    if (VERBOSE_MODE >= 1)
      fprintf('\n%s(K=%d) optimizing trial %d: %g', Verbose_Prefix, K, my_it, my_LL);
    end
    
    % estimate hyperparameters
    
    % copy options
    % v0.80 ABC - remove trial initialization field to save space
    vbhemopt2 = rmfield(vbhemopt, 'VHintial');
    vbhemopt2.initmode = 'inith3m';
    vbhemopt2.inithmm  = my_h3m;
    
    %       h3m_news{my_it} = vbhem_h3m_c_hyp(h3m_b, vbhemopt2);
    %       LLall(my_it)  = h3m_news{my_it}.LL;
    
    % run vbhem co-clustering
    tmpout = vbhem_coh3m_c_hyp(h3m_bs, vbhemopt2);
    
    % v0.80 ABC - cleanup output
    % ...remove initial h3m to save space
    tmpout.learn_hyps.inithmm = 'removed';
    tmpout.learn_hyps.vbopt.inithmm = 'removed';
    
    h3m_news1{q} = tmpout;
    LLall_1(q)  = h3m_news1{q}.LL;
    
    % visualization
    if (show_LL_plot)
      figure(VLL)
      legs4 = plot(my_it, LLall(my_it), 'rx');
      legend([legs, legs4], {'original LL', 'selected trial', 'optimized LL', 'ignore bounds'});
      drawnow
    end
    
  end
  
  
  h3m_news = cell(sizeh3m_news);
  h3m_news(unique_h3m_it) = h3m_news1(:,1);
  LLall(unique_h3m_it) = LLall_1;
  clear h3m_news1
  
end % end do_learn_hyps


%% choose the best
[maxLL,maxind] = max(LLall);
h3m_new = h3m_news{maxind};
% if (vbhemopt.keep_ini_hmms)
%     h3m_new.inih3m = All_inihmms{maxind};
% end
t_best = maxind;
LL_best = h3m_new.LL;  


if (VERBOSE_MODE >= 1)  
    fprintf('\nBest run is %d: LL=%g\n\n',t_best,LL_best)
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
    h3m_new.suboptimal_hmms = h3m_news; %unique_h3m;
end





