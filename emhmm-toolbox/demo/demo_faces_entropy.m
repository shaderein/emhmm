% demo_faces_entropy - example computing the entropy of HMMs
% 
% assumes that demo_faces has already been run
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-07-21
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% 2021-07-21: v0.78 - initial version

clear
close all

%% Data files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% names of files used
xlsname = 'demodata.xls';    % Excel File with fixation data
faceimg = 'ave_face120.png';  % average face image

% names of files for saving results
matfile_individual = 'models_demo_faces_individual.mat';
matfile_group      = 'models_demo_faces_group.mat';



%% load individual models from disk
fprintf('loading individual models to %s\n', matfile_individual);
load(matfile_individual)

% load data
load('demodata.mat')

%% for demo, only look at 3 examples
ss = [1,4,5];
hmms = hmms(ss);
data = data(ss);


%% compute entropy values
T = 5;  % maximum length to compute

% for each HMM...
for i=1:length(hmms)
  % select the ith HMM
  myhmm = hmms{i};
  
  % overall entropy
  ent(i) = vbhmm_entropy(myhmm, data{i});
  
  % sequence joint entropy from 1:tt (increasing lengths)
  for tt=1:T
    entj(i,tt) = vbhmm_entropy_joint(myhmm, tt);
  end
  
  % marginal entropy at time t
  entm(i,:) = vbhmm_entropy_marginal(myhmm, T);
  
  % conditional entropy at time t given t-1
  entc(i,1) = nan; % can't compute for first fixation
  for tt=2:T
    entc(i,tt) = vbhmm_entropy_conditional(myhmm, tt);
  end
  
  % conditional entropy (steady-state) t=inf
  entc_ss(i) = vbhmm_entropy_conditional(myhmm);
  
  % joint entropy of fixation pairs (t-1,t)
  entj2(i,1) = nan; % can't compute for first fixation
  for tt=2:T
    entj2(i,tt) = vbhmm_entropy_joint2(myhmm, tt);
  end
  
  %  mutual information of fixation pairs (t-1,t)
  mi(i,1) = nan; % can't compute for first fixation
  for tt=2:T
    mi(i,tt) = vbhmm_mutualinfo2(myhmm, tt);
  end
end



%% visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup sampling step for density heat maps
img = imread(faceimg);
imgy = size(img,1);
imgx = size(img,2);
Xp = 1:imgx;
Yp = 1:imgy;

figure(1)

% for each HMM
for i=1:length(hmms)
  myhmm = hmms{i};
  
  
  % plot HMM
  ax=subplot(length(hmms),5,5*(i-1)+1);
  vbhmm_plot_compact(myhmm, img);
  title(sprintf('HMM %d\noverall ent=%0.3g', i, ent(i)));  

  % compute marginal probability distributions
  pym = vbhmm_prob_marginalobs(myhmm, T);
  pxm = vbhmm_prob_marginalstate(myhmm, T);
    
  % plot each marginal distribution as a heatmap
  for t=1:4
    pym_h{t} = gmm2heatmap(pym{t}, Xp, Yp);    
    ax = subplot(length(hmms),5,5*(i-1)+t+1);
    imagesc(Xp, Yp, pym_h{t})
    colormap(ax, jet)
    axis image
    axis off
    title(sprintf('t=%d\nmarg ent=%0.3g', t, entm(i,t)))
  end
    
end

% plot graphs of entropy vs time
figure(2)
subplot(2,2,1)
plot(1:T, entm, 'x-');
xlabel('t');
ylabel('H');
grid on
title('marginal entropy')
legs = {};
for i=1:length(hmms)
  legs{i} = sprintf('HMM %d', i);
end
legend(legs)


subplot(2,2,2)
plot(2:T, entc(:,2:end), 'o-');
xlabel('t');
ylabel('H');
xlim([1,5])
grid on
title('conditional entropy')

subplot(2,2,3)
plot(2:T, entj2(:,2:end), 'o-');
xlabel('t');
ylabel('H');
title('pairwise joint entropy');
xlim([1,5])
grid on

subplot(2,2,4)
plot(2:T, mi(:,2:end), 'o-');
xlabel('t');
ylabel('I');
title('pairwise mutual information');
xlim([1,5])
grid on
