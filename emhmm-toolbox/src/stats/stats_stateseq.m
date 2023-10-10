function [sort_seqs, sort_pr] = stats_stateseq(hmm, seqlen)
% stats_stateseq - calculate the highest probability hidden-state (ROI) sequences
%
%  [seqs, pr] = stats_stateseq(hmm, seqlen)
%
% INPUTS
%     hmm = HMM learned with vbhmm_learn or vbhem_cluster
%  seqlen = ROI sequence length to use (default=3)
%
% OUTPUTS
%    seqs = the ROI sequences, sorted by descending probability
%      pr = the probabilities for each sequence
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-05-25
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

% v0.73 - rename vbhmm_state_prob to vbhmm_prob_state

if nargin<2
  seqlen = 3;
end

K = length(hmm.pdf);

% get all possible state sequences
tmp = cell(1,seqlen);
[tmp{:}] = ndgrid(1:K);
for i=1:seqlen
  tmp{i} = tmp{i}(:);
end
tmp2 = cat(2,tmp{:});
seqs = mat2cell(tmp2, ones(size(tmp2,1),1), seqlen);
N = length(seqs);

% calculate probabilites
pr = vbhmm_prob_state(hmm, seqs);

% sort by probability
[sort_pr, sort_i] = sort(pr(:), 1, 'descend');
sort_seqs = {seqs{sort_i}};


