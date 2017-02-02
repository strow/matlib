clear dbt dbt_err
for i=1:40
   for j=2:8
      fn = ['fit' int2str(j) '_robs_lat' int2str(i)];
      g = load(fn);
      dbt(i,j,:) = g.dbt;
      dbt_err(i,j,:) = g.dbt_err;
   end
end