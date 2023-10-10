function err = totalerror(Ain, Bin, mode, verbose, indentsp)
% totalerror - compute the total error between A and B
%
%    err = totalerror(A, B, mode, verbose)
%
%        mode = 1 (abs error)  default
%             = 2 (squared error)
%             = 3 (max abs error)
%             = 4 (max percent change, A=new, B=old)
%     verbose = 0 - say nothing (default)
%               1 - list differences
%               2 - list differences when not zero
%
%    A and B can be matrices, cells, structs, or any combo thereof
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-08-08
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

if (nargin<3)
  mode = 1;
end
if (nargin<4)
  verbose = 0;
end
if (nargin<5)
  indentsp = 0;
end

myindent = repmat('  ', [1 indentsp]);

if (mode == 3 || mode == 4)
  maxmode = 1;
else
  maxmode = 0;
end

% recurse on cell contents
if (iscell(Ain) && iscell(Bin))
  if (~all(size(Ain)==size(Bin)))
    error('totalerror: cells not the same size!');
  end
  err = 0;
  for i=1:numel(Ain)
    if (verbose)
      fprintf('%scell{%d}:\n', myindent, i);
    end
    myerr = totalerror(Ain{i}, Bin{i}, mode, verbose, indentsp+1);
    if (maxmode)
      err = max(err, myerr);
    else
      err = err + myerr;
    end
    if (verbose == 1) || ((verbose == 2) && (myerr ~= 0))
      fprintf('%scell{%d}: totalerror = %g\n', myindent, i, myerr);
    end
  end
  
% recurse on struct contents
elseif (isstruct(Ain) && (isstruct(Bin)))
  fnames = fieldnames(Ain);
  err = 0;
  for i=1:length(fnames)
    myname = fnames{i};
    if (isfield(Bin, myname))
      if (verbose)
	fprintf('%sstruct "%s":\n', myindent, myname);
      end
      myerr = totalerror(getfield(Ain, myname), ...
			 getfield(Bin, myname), ...
			 mode, verbose, indentsp+1);
      if (verbose == 1) || ((verbose == 2) && (myerr ~= 0))
	fprintf('%sstruct "%s": totalerror = %g\n', myindent, myname, myerr);
      end
      if (maxmode)
	err = max(err, myerr);
      else
	err = err + myerr;
      end
    else
      warning(sprintf('totalerror: B does not have field "%s"', myname));
    end
  end
  
  fnames = fieldnames(Bin);
  for i=1:length(fnames)
    myname = fnames{i};
    if (~isfield(Ain, myname))
      warning(sprintf('totalerror:A does not have field "%s"', myname));
    end
  end
  
% compare 2 matrices
else
  
  Asize = size(Ain);
  Bsize = size(Bin);
  if (~all(Asize == Bsize))
    warning('matrices are different sizes, truncating');
    newsize = min(Asize, Bsize);
    switch( length(newsize))
     case 2
      Ain = Ain(1:newsize(1), 1:newsize(2));
      Bin = Bin(1:newsize(1), 1:newsize(2));
     case 3
      Ain = Ain(1:newsize(1), 1:newsize(2), 1:newsize(3));
      Bin = Bin(1:newsize(1), 1:newsize(2), 1:newsize(3));
     case 4
      Ain = Ain(1:newsize(1), 1:newsize(2), 1:newsize(3), 1:newsize(4));
      Bin = Bin(1:newsize(1), 1:newsize(2), 1:newsize(3), 1:newsize(4));      
     otherwise
      warning('ndims > 4 not supported!');
    end
  end
  
  switch(mode)
   case 1
    D = abs(Ain-Bin);
    err = sum(D(:));
  
   case 2
    D = (Ain-Bin).^2;
    err = sum(D(:));
    
   case 3
    D = abs(Ain-Bin);
    err = max(D(:));
    
   case 4
    D = abs(1-Bin./Ain);
    err = max(D(:));
    
   otherwise
    error('unknown mode');
  end
   


end

  
