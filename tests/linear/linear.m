
% construct a list of seeds.
seeds = randi([11111, 99999], 200, 1);

% initialize the output file.
fh = fopen('linear.dat', 'w');
D = [];

% loop over the seeds.
for n = 1 : length(seeds),
  % build the command string.
  cmd = ['env DYLD_LIBRARY_PATH=../../lib', ...
            ' RNG_SEED=', num2str(seeds(n)), ...
            ' ./linear | tail -n 5'];

  % run the command.
  [status, output] = system(cmd);
  dat = str2num(output);

  % get the frequencies.
  f1 = dat(:,1);
  f2 = dat(:,3);

  % get the weights.
  w1 = dat(:,2);
  w2 = dat(:,4);

  % get the extras.
  T = dat(:,5:end);

  % loop over the signals.
  while (any(!isnan(f1)))
    % find the current best match.
    df = (f1 - f2').^2;
    [i, j] = find(df == min(vec(df)));

    % output the matching values in the result.
    fprintf(fh, '%16.9e %16.9e %16.9e %16.9e', ...
            f1(i), f2(j), w1(i), w2(j));

    % output the extras.
    for tc = 1 : columns(T)
      fprintf(fh, ' %16.9e', T(j,tc));
    end
    fprintf(fh, '\n');

    % avoid making repeat matches.
    f1(i) = nan;
    f2(j) = nan;
  end
  fflush(fh);
end

% close the output file.
fclose(fh);

