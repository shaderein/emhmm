function h3m_r_new = hem_fix_degenerate_component(h3m_r_new, i_z, N, Nr)
% hem_fix_degenerate_component - fix a degenerate component [internal function]
%
% h3m_r_new = hem_fix_degenerate_component(h3m_r_new, i_z, N, Nr)
%
% i_z is the index of the component to fix
%
% ---
% H3M Toolbox 
% Copyright 2012, Emanuele Coviello, CALab, UCSD
% Antoni Chan 2017/06/06

%fprintf('!!! modifying h3m: one hmm has zero prior \n')
[foo highest] = max(h3m_r_new.omega);
h3m_r_new.omega([i_z highest]) = h3m_r_new.omega(highest)/2;

% renormalize for safety
h3m_r_new.omega = h3m_r_new.omega / sum(h3m_r_new.omega);
h3m_r_new.hmm{i_z} = h3m_r_new.hmm{highest};

% perturb
h3m_r_new.hmm{i_z}.prior = h3m_r_new.hmm{highest}.prior + (.1/N) * rand(size(h3m_r_new.hmm{highest}.prior));
A = h3m_r_new.hmm{highest}.A;
f_zeros = find(A == 0);
A = (.1/N) * rand(size(A));
A(f_zeros) = 0;
        
h3m_r_new.hmm{i_z}.A     = A;
% renormalize
h3m_r_new.hmm{i_z}.prior = h3m_r_new.hmm{i_z}.prior / sum(h3m_r_new.hmm{i_z}.prior);
%% BUG FIX: 2017-05-25: ABC: N -> Nr
h3m_r_new.hmm{i_z}.A     = h3m_r_new.hmm{i_z}.A    ./   repmat(sum(h3m_r_new.hmm{i_z}.A,2),1,Nr);

