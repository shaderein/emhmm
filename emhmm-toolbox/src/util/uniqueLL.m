function [inds, VLL, legs] = uniqueLL(LLall, diffthresh, show_LL_plot)
% uniqueLL - get indicies of unique LL within a threshold
%
%  [inds, fignum, legs] = uniqueLL(LLs, diffthresh, show_LL_plot)
%
%  INPUTS
%            LLs - list of log-likelihoods
%     diffthresh - difference threshold for LL to be similar
%   show_LL_plot - 1 = show a plot [default=0]
%  
% OUTPUTS
%      inds = indices for unique LLs
%    fignum = handle to figure created with show_LL_plot
%      legs = legend handle for figure
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-01-13
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong


legs = [];
VLL = [];
numits = length(LLall);

unique_hmm_LL = [];
unique_hmm_it = [];
              
% visualization of LLs
if (show_LL_plot)
  VLL = figure;
  legs1 = plot(1:length(LLall), LLall, 'x');
  xlabel('trial');
  ylabel('LL');
  hold on
  grid on
end

% go through the HMMs, select unique ones for hyp estimation
for it=1:numits
  my_LL = LLall(it);
  
  new_hmm = 0;
  if (it==1)
    new_hmm = 1;
  else
    % check if LL is similar to previous models
    pLL =abs((unique_hmm_LL - my_LL) / my_LL);
    if (abs(pLL) > diffthresh)
      new_hmm = 1;
      
    else
      % it's possible that two HMMs will have the same LL but different configurations.
      % ...we could check if the parameters are unqiue enough.
    end
  end
  
  % found a unique HMM
  if (new_hmm)
    % add to the list
    unique_hmm_LL(end+1) = my_LL;
    unique_hmm_it(end+1) = it;
                        
    if (show_LL_plot)
      figure(VLL)
      tmphi = my_LL / (1-diffthresh);
      tmplo = my_LL / (1+diffthresh);      
      %legs3 = plot([1 numits], [tmplo tmplo], 'b-');
      %plot([1 numits], [tmphi tmphi], 'b-');
      legs3 = fill([1 1 numits numits 1], [tmplo tmphi tmphi tmplo tmplo], [0.9 0.9 0.9], ...
        'facealpha', 0.5, 'edgecolor', 'none');
      legs2 = plot(it, my_LL, 'ro');
      legs = [legs1, legs2, legs3];
      legend(legs, {'original LL', 'selected trial', 'ignore bounds'});
      drawnow
    end
            
  end
end

inds = unique_hmm_it;