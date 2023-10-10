function [hmm,L] = vbhmm_em_hyp(data, K, vbopt)
% vbhmm_em_hyp - estimate hyperparameters for vbhmm (internal function)
% Use vbhmm_learn instead. For options, see vbhmm_learn.
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-07-31
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% uses vbopt.learn_hyps option
% requires initmode='inithmm', and inithmm is specified

% see vbhmm_learn for changes

%% set default hyps
if ~iscell(vbopt.learn_hyps)
  if (vbopt.learn_hyps == 1)
    learn_hyps = {'alpha0', 'epsilon0', 'v0', 'beta0', 'W0', 'mu0'};
  else
    error('learn_hyps not set properly');
  end
else
  learn_hyps = vbopt.learn_hyps;
end

%% select the hyps for optimization
hypinfo = get_hypinfo(learn_hyps, vbopt);

% check if we already have an initial HMM
if ~strcmp(vbopt.initmode, 'inithmm')
  error('requires initization using inithmm');
end

%% initial parameter vector, and the objective function
initX = init_hyp(vbopt, hypinfo);
myf   = @(X) vbhmm_grad(data, K, set_vbopt(X, vbopt, hypinfo), hypinfo);

%vbopt.verbose = 2;
myoutfun = @(x,optimValues,state) outfun(x,optimValues,state,vbopt.verbose);

switch(vbopt.minimizer)
  case 'fminunc'
    %% do the optimization using fminunc

    % set optimization options
    options = optimoptions(...
      'fminunc', ...
      'Algorithm', 'trust-region', ...
      'GradObj',   'on', ...
      'Display', 'off', ...
      'TolFun',    vbopt.minDiff);
    %'Display',   'Iter-Detailed',
        
    if (vbopt.verbose >= 1)
      options = optimoptions(options, 'OutputFcn', myoutfun);
    end
    
    % optimize
    [opt_transhyp, opt_L, exitflag] = fminunc(myf, initX, options);
    
    %opt_transhyp'
    %opt_L
    
    % check for problem in exitflag
    if (exitflag==0)
      warning('not enough iterations to fully optimize');
    elseif (exitflag<0)
      error('error optimizing the function');
    end

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
    
    %fprintf('%g', opt_L);    
    %[opt_transhyp, opt_Ls] = minimize(initX, myf, 100);  
    
  otherwise
    error('bad minimizer specified');
end

%% run EM to get the final model again
% (since fminunc doesn't return results from function evaluations)
vbopt2 = set_vbopt(opt_transhyp, vbopt, hypinfo);
%hmm = vbhmm_learn(data, K, vbopt2);
hmm = vbhmm_em(data, K, vbopt2);
L = hmm.LL;

% dump optimized hyps
if (vbopt.verbose>1)
  fprintf('\n  ');
  vbhmm_print_hyps({hypinfo.optname}, vbopt2);
end

%% add info about learning hyps
hmm.learn_hyps.hypinfo      = hypinfo;
hmm.learn_hyps.opt_transhyp = opt_transhyp;
hmm.learn_hyps.opt_L        = -opt_L;  % opt_L is actually the negative Log-likelihood

% v0.74 - remove learn_hyps.inithmm - it's the same as learn_hyps.vbopt.inithmm
%hmm.learn_hyps.inithmm      = vbopt2.inithmm;
hmm.learn_hyps.vbopt        = vbopt2;

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

%% function for calculating the objective and gradient %%%%%%%%%%%%%%%%%%%
function [L, dL] = vbhmm_grad(data, K, vbopt, hypinfo)
% X is the parameter vector

% calculate derivative?
if (nargout>1)
  vbopt.calc_LLderiv = 1;
end

% run EM
hmm = vbhmm_em(data, K, vbopt);
L  = -hmm.LL;

if (nargout>1)  
  dL = [];
  % append gradients
  for i=1:length(hypinfo)
    dL = [dL; -hmm.dLL.(hypinfo(i).derivname)(:)];
  end
end


if (vbopt.verbose >=2)
  fprintf('** gradient: dL=');
  fprintf('(');
  fprintf('%g ', -dL);
  fprintf(')\n');
end

% DEBUGGING: check for invalid gradients
%if any(isnan(dL))
%  keyboard
%end
%if any(~isreal(dL))
%  keyboard
%end
  

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


