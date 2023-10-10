function [roiseq] = vbhmm_map_state(hmm, data)
% vbhmm_map_state - compute the most probable state (ROI) sequence from fixation data
%
%  [roiseq] = vbhmm_map_state(hmm, data)
%
% INPUTS
%   hmm = HMM learned with vbhmm_learn
%  data = 1xN cell array of fixation sequences (as in data passed to vbhmm_learn)
%
% OUTPUTS
%   roiseq = 1xN cell array of corresponding ROI seequences
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-08-22
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

if ~iscell(data)
  data = {data};
end

K = length(hmm.prior);
N = length(data);

roiseq = cell(size(data));

for j = 1:N
  temp = data{j};
  emit_p = zeros(size(temp,1), K);
  for k = 1:K
    emit_p(:,k) = mvnpdf(temp, hmm.pdf{k}.mean, hmm.pdf{k}.cov);
  end
  [path,loglik] = viterbi_path(hmm.prior, hmm.trans, emit_p');
  roiseq{j} = path;
end

end


function [path,loglik] = viterbi_path(prior, transmat, obslik)
% VITERBI Find the most-probable (Viterbi) path through the HMM state trellis.
% path = viterbi(prior, transmat, obslik)
%
% Inputs:
% prior(i) = Pr(Q(1) = i)
% transmat(i,j) = Pr(Q(t+1)=j | Q(t)=i)
% obslik(i,t) = Pr(y(t) | Q(t)=i)
%
% Outputs:
% path(t) = q(t), where q1 ... qT is the argmax of the above expression.


% delta(j,t) = prob. of the best sequence of length t-1 and then going to state j, and O(1:t)
% psi(j,t) = the best predecessor state, given that we ended up in state j at t

scaled = 1;

T = size(obslik, 2);
prior = prior(:);
Q = length(prior);

delta = zeros(Q,T);
psi = zeros(Q,T);
path = zeros(1,T);
scale = ones(1,T);


t=1;
delta(:,t) = prior .* obslik(:,t);
if scaled
  [delta(:,t), n] = normalize(delta(:,t));
  scale(t) = 1/n;
end
psi(:,t) = 0; % arbitrary value, since there is no predecessor to t=1
for t=2:T
  for j=1:Q
    [delta(j,t), psi(j,t)] = max(delta(:,t-1) .* transmat(:,j));
    delta(j,t) = delta(j,t) * obslik(j,t);
  end
  if scaled
    [delta(:,t), n] = normalize(delta(:,t));
    scale(t) = 1/n;
  end
end
[p, path(T)] = max(delta(:,T));
for t=T-1:-1:1
  path(t) = psi(path(t+1),t+1);
end

% If scaled==0, p = prob_path(best_path)
% If scaled==1, p = Pr(replace sum with max and proceed as in the scaled forwards algo)
% Both are different from p(data) as computed using the sum-product (forwards) algorithm

%if 0
if scaled
  loglik = -sum(log(scale));
  %loglik = prob_path(prior, transmat, obslik, path);
else
  loglik = log(p);
end
%end
end


function [M, z] = normalize(A, dim)
% NORMALISE Make the entries of a (multidimensional) array sum to 1
% [M, c] = normalise(A)
% c is the normalizing constant
%
% [M, c] = normalise(A, dim)
% If dim is specified, we normalise the specified dimension only,
% otherwise we normalise the whole array.

if nargin < 2
  z = sum(A(:));
  % Set any zeros to one before dividing
  % This is valid, since c=0 => all i. A(i)=0 => the answer should be 0/1=0
  s = z + (z==0);
  M = A / s;
elseif dim==1 % normalize each column
  z = sum(A);
  s = z + (z==0);
  %M = A ./ (d'*ones(1,size(A,1)))';
  M = A ./ repmatC(s, size(A,1), 1);
else
  % Keith Battocchi - v. slow because of repmat
  z=sum(A,dim);
  s = z + (z==0);
  L=size(A,dim);
  d=length(size(A));
  v=ones(d,1);
  v(dim)=L;
  %c=repmat(s,v);
  c=repmat(s,v');
  M=A./c;
end
end

