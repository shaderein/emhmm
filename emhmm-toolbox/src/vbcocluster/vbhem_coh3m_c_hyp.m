function [h3m,L] = vbhem_coh3m_c_hyp(h3m_bs, vbopt)
% vbhem_h3m_c_hyp - estimate hyperparameters for vbh3m_cluster (internal function)
%  For options, see vbh3m_cluster.
%
% uses vbopt.learn_hyps option
% requires initmode='inith3m', and inith3m is specified
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-07-29
% Antoni B. Chan, Janet H. Hsiao, Lan Hui
% City University of Hong Kong, University of Hong Kong

% 2021-07-29: v0.80 - Hui - initial version


%% set default hyps
if ~iscell(vbopt.learn_hyps)
  if (vbopt.learn_hyps == 1)
    learn_hyps = {'alpha0','eta0' 'epsilon0', 'v0', 'lambda0', 'W0', 'm0'};
  else
    error('learn_hyps not set properly');
  end
else
  learn_hyps = vbopt.learn_hyps;
end

%% select the hyps for optimization
hypinfo = vbhem_get_hypinfo(learn_hyps, vbopt);

% check if we already have an initial H3M
if ~strcmp(vbopt.initmode, 'inith3m')
  error('requires initization using inith3m');
end

%% initial parameter vector, and the objective function
initX = init_hyp(vbopt, hypinfo);
myf   = @(X) vbh3m_cograd(h3m_bs, set_vbopt(X, vbopt, hypinfo), hypinfo);

myoutfun = @(x,optimValues,state) outfun(x,optimValues,state,vbopt.verbose);

switch(vbopt.minimizer)

  case {'minimize-lbfgs', 'minimize-bfgs', 'minimize-cg'}
    %% do the optimization using minimize_new (LBFGS)
    p.length = 100;
    if strcmp(vbopt.minimizer, 'minimize-lbfgs')
      p.method = 'LBFGS';
    elseif strcmp(vbopt.minimizer, 'minimize-cg')
      p.method = 'CG';
    else
      p.method = 'BFGS';
    end    
    p.verbosity = vbopt.verbose;
    p.outfun    = myoutfun;
    [opt_transhyp, opt_Ls] = minimize_new(initX, myf, p);
    opt_L = opt_Ls(end);
    
  otherwise
    error('bad minimizer specified');
end

%% run EM to get the final model again
vbopt3 = set_vbopt(opt_transhyp, vbopt, hypinfo);
%hmm = vbhmm_learn(data, K, vbopt2);
[vbopt2, hyp_clipped] = vbhem_clip_hyps(vbopt3);
vbopt2.hyp_clipped = hyp_clipped;

%initialize all stimuli
h3ms =[];
moptpointer = vbopt2;
for i=1:vbopt2.stimulisize
    if length(vbopt2.S)>1
        moptpointer.S = vbopt2.S(i);
    end
    moptpointer.inithmm = vbopt2.inithmm.h3m_rs{i};
    h3m1 = vbhemhmm_init(h3m_bs(i),moptpointer);
    h3ms = [h3ms; h3m1];
end


h3m = vbhem_coh3m_c_step_fc(h3ms,h3m_bs,vbopt2);

L = h3m.LL;
h3m.iter = vbopt2.inithmm.iter; 
% dump optimized hyps
if (vbopt.verbose>1)
  fprintf('\n  ');
  vbhmm_print_hyps({hypinfo.optname}, vbopt2);
end

%% add info about learning hyps
h3m.learn_hyps.hypinfo      = hypinfo;
h3m.learn_hyps.opt_transhyp = opt_transhyp;
h3m.learn_hyps.opt_L        = -opt_L;  % opt_L is actually the negative Log-likelihood
h3m.learn_hyps.inithmm      = vbopt2.inithmm;
h3m.learn_hyps.vbopt        = vbopt2;

% show subject
if 0
  vbhmm_plot(opt_hmm, mydata, faceimg);
  figure,
  subplot(2,1,1)
  vbhmm_plot_compact(hmm, faceimg);
  title(sprintf('original LL=%g', hmm.LL));
  subplot(2,1,2)
  vbhmm_plot_compact(opt_hmm, faceimg);
  title(sprintf('hyp optimized LL=%g', opt_hmm.LL));

  %% plot new m
  if any(strcmp(optname, 'mu'))
    figure
    img = imread(faceimg);
    plot_emissions(mydata, opt_hmm.gamma, opt_hmm.pdf, img)
    hold on
    legs(1) = plot(vbopt3.mu(1), vbopt3.mu(2), 'kx', 'Markersize', 10);
    legs(2) = plot(vbopt.mu(1), vbopt3.mu(2), 'ko', 'Markersize', 10);
    meanfix = mean(cat(1,mydata{:}));
    legs(3) = plot(meanfix(1), meanfix(2), 'k+', 'Markersize', 10);
    hold off
    legend(legs, {'optimized', 'old', 'mean'});
  end
end


%% function to modify vbopt using parameter vector X %%%%%%%%%%%%%%%%%%%%%%
function vbopt2 = set_vbopt(X, vbopt, hypinfo)
vbopt2 = vbopt;
ind = 1;
for i=1:length(hypinfo)
  myinds = ind:(ind+hypinfo(i).hypdims-1);
  transhyp = hypinfo(i).hyptrans(X(myinds));  % transform from opt-space to hyp-space
  vbopt2.(hypinfo(i).optname) = transhyp;
  ind = ind+hypinfo(i).hypdims;
  
  %if any(~isreal(transhyp))
  %  keyboard
  %end

end

if (vbopt.verbose >=2)
  fprintf('Setting vbopt: ');
  fprintf('X=(');
  fprintf('%g ', X);
  fprintf(')\n');
end
  

%% function for making the initial parameter vector %%%%%%%%%%%%%%%%%%%%%%%%
function [X] = init_hyp(vbopt, hypinfo)
X = [];
for i=1:length(hypinfo)
  newX = hypinfo(i).hypinvtrans( vbopt.(hypinfo(i).optname) );  
  X = [X; newX(:)];
end

%if any(~isreal(X))
%  keyboard
%end


%% function for displaying progress
function stop = outfun(x,optimValues,state, verbose)
persistent lastfval;

stop = false;
switch state
  case 'init'
    % do nothing
    lastfval = optimValues.fval;
  case 'iter'
    if (verbose == 1)
      fprintf('.');
    elseif (verbose > 1)
      % debugging stuff
      dLL = (lastfval-optimValues.fval) / lastfval;
      fprintf('** optimization iteration **\n');
      fprintf('LL=%g [dLL=%g] ', -optimValues.fval, dLL);
      fprintf('(');
      fprintf('%g ', x);
      fprintf(')\n');
    end
    lastfval = optimValues.fval;
    
  case 'done'
    if (verbose==1)
      fprintf('LL=%g', -optimValues.fval);
    elseif (verbose>1)
      fprintf('** optimized LL=%g\n', -optimValues.fval);
    end
  otherwise
    % do nothing
end


