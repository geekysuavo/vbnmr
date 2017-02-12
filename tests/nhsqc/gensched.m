
n = 300;

rates = [0.005, 0.02];

Tmin = [0, 0];
Tmax = [0.02, 0.1];

dat = [];

while rows(dat) < n
  k = randperm(4, 1) - 1;
  s = arrayfun(@(x) exprnd(x), rates);

  if all(s >= Tmin) && all(s <= Tmax)
    dat = [dat; k, s, 0];
  end
end

fh = fopen('sched.dat', 'w');

for ii = 1 : rows(dat)
  fprintf(fh, '%d %e %e %e\n', ...
          dat(ii,1), dat(ii,2), ...
          dat(ii,3), dat(ii,4));
end

fclose(fh);

