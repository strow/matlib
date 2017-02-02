function h = pbias(f,bias,bias_std,count,lati,dayi)
   
count = squeeze(count(dayi,lati,:));
count = count./max(count);
mbias = squeeze(bias(dayi,lati,:));
mstd  = squeeze(bias_std(dayi,lati,:));

ig = find(count > 0.98 * max(count) &  abs(mstd) < 5  );

h = plot(f(ig),mbias(ig));
hold on;
plot(f(ig),mstd(ig))

