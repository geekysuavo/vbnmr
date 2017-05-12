
% set the number of signals per experiment
% (must match linear.c / linear.vfl)
M = 5;

% open the input file.
dats = load(argv(){1});

% loop over the experiments.
for idx = 1 : M : length(dats)
  % get the experiment data.
  dat = dats(idx : idx + M - 1, :);

  % get the initial parameters.
  f1 = dat(:,1);
  w1 = dat(:,2);

  % get the inferred parameters.
  f2 = dat(:,3);
  w2 = dat(:,4);

  % get the extras.
  T = dat(:,5:end);

  % loop over the signals.
  while (any(!isnan(f1)))
    % find the current best match.
    df = (f1 - f2').^2;
    [i, j] = find(df == min(vec(df)));

    % output the matching values in the result.
    printf('%16.9e %16.9e %16.9e %16.9e', ...
           f1(i), w1(i), f2(j), w2(j));

    % output the extras.
    for tc = 1 : columns(T)
      printf(' %16.9e', T(j,tc));
    end
    printf('\n');

    % avoid making repeat matches.
    f1(i) = nan;
    f2(j) = nan;
  end
end

