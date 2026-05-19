function printarray(raaX,comment,iN)

if nargin == 1
  iN = 10;
  comment = [];
elseif nargin == 2
  iN = 10;
end

[mm,nn] = size(raaX);
if mm == 1 | nn == 1
  raaX = reshape(raaX,1,length(raaX));
  [mm,nn] = size(raaX);
end

[mm,nn] = size(raaX);
len = length(raaX);
len = nn;

for kk = 1 : mm
  raX = raaX(kk,:);  
  if isa(raX,'integer') | sum(abs(raX - floor(raX))) == 0 
    iaIorR(kk) = +1;
  elseif isa(raX,'single') | isa(raX,'double') 
    iaIorR(kk) = -1;
  end
end
if mm == 1
  iIorR = iaIorR;
else
  iIorR = -1;
  if sum(iaIorR) == mm
    iIorR = +1;
  end
end


if length(comment) > 0
  if  nn <= iN & mm == 1
    fprintf(1,'%s \n',comment)
  else
    fprintf(1,'input matrix of size %10i x %10i : %s \n',mm,nn,comment)
  end
end

if nn <= iN
  for kk = 1 : mm
    raX = raaX(kk,:);
    if iIorR == +1
      for ii = 1 : length(raX)
        fprintf(1,'%8i ',raX(ii));
      end
    elseif iIorR == -1
      for ii = 1 : length(raX)
        fprintf(1,'%12.6f ',raX(ii));
      end
    end
    fprintf(1,' \n');
  end
else
  iaN = (1:iN);
  iRepeat = floor(length(raaX)/iN);
  for kk = 1 : iRepeat
    iaInd = iaN + (kk-1)*iN;
    for jj = 1 : mm
      %raX = raaX(iaInd);   
      raX = raaX(jj,iaInd);
      if iIorR == +1
        for ii = 1 : iN
          fprintf(1,'%8i ',raX(ii));
        end
      elseif iIorR == -1
        for ii = 1 : iN
          fprintf(1,'%12.6f ',raX(ii));
        end
      end
      fprintf(1,' \n');
    end
    if mm > 1
      fprintf(1,' \n');
    end
  end

  jj = iRepeat + 1;
  iaInd = iaN + (jj-1)*iN;
  boo = find(iaInd <= length(raaX));
  iaInd = iaInd(boo);
  for jj = 1 : mm
    %raX = raaX(iaInd);   
    raX = raaX(jj,iaInd);
    if iIorR == +1
      for ii = 1 : length(raX)
        fprintf(1,'%8i ',raX(ii));
      end
    elseif iIorR == -1
      for ii = 1 : length(raX)
        fprintf(1,'%12.6f ',raX(ii));
      end
    end
    fprintf(1,'\n')
  end
  if mm > 1
    fprintf(1,'\n')
  end
end

