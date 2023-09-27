function hmmi = vbhmm_group2ind(hmmg)
% vbhmm_group2ind - convert a HMM group into individual HMMs
% 

G = length(hmmg.group_ids);

hmmi = cell(1,G);
for g=1:G
  myinds = hmmg.group_inds{g};
  
  % copy all first
  hmmi{g} = hmmg;
  
  % recopy some things.
  hmmi{g}.prior  = hmmg.prior{g};
  hmmi{g}.trans  = hmmg.trans{g};
  hmmi{g}.gamma  = hmmg.gamma(myinds);
  hmmi{g}.M      = hmmg.M{g};
  hmmi{g}.N1     = hmmg.N1{g};
  hmmi{g}.N      = hmmg.Ng{g};
  hmmi{g}.mygroup_ids  = hmmg.group_ids(g);
  hmmi{g}.mygroup_inds = hmmg.group_inds{g};
  hmmi{g}.varpar.epsilon = hmmg.varpar.epsilon{g};
  hmmi{g}.varpar.alpha   = hmmg.varpar.alpha{g};
  
  % remove things
  hmmi{g} = rmfield(hmmi{g}, {'Ng', 'group_map', 'group_ids', 'group_inds'});
end

