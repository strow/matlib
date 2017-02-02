btobsx = btobs;

for i=1:352
   for j=1:40
      ig = goodchan_stats(squeeze(count(i,j,:)));
      ib = setxor(ig,1:2378);
      btobsx(i,j,ib) = NaN;
   end
end
