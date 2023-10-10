function mix = naiveWordMix(allMix)
%allMix is an array of gaussian mixture
%mix is the resulting mixture after naive model averaging

% 2017-08-16: ABC: full cov support

numMixs = size(allMix,2);

ncentres = 0;
priors = [];
centres = [];
covars = [];
nwts = 0;

for i = 1:numMixs
  ncentres     = ncentres + allMix(i).mix.ncentres;
  priors       = [priors   allMix(i).mix.priors];
  centres      = [centres;  allMix(i).mix.centres];
  switch(allMix(1).mix.covar_type)
    case 'diag'
      covars       = [covars;   allMix(i).mix.covars];
    case 'full'
      covars       = cat(3, covars, allMix(i).mix.covars);
  end
  if ~isfield(allMix(i).mix,'nwts')
    nwts         = nwts     + 1;
  else
    nwts         = nwts     + allMix(i).mix.nwts;
  end
end

mix.type = allMix(1).mix.type;
mix.nin = allMix(1).mix.nin;
mix.ncentres = ncentres;
mix.covar_type = allMix(1).mix.covar_type;
mix.priors = priors ./ repmat(numMixs,1, ncentres);
mix.centres = centres;
mix.covars = covars;
mix.nwts = nwts;