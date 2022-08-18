function dump_result(data, path, prefix)
% Get a pairwise distance stats
pos = cell2mat(arrayfun(@(x) x.pos, data, 'UniformOutput', false));
sur = arrayfun(@(x) x.surface, data);
vol = arrayfun(@(x) x.volume,  data);
rat = arrayfun(@(x) x.ratio,   data);

cmass_dist = pdist2(pos, mean(pos));

pairwise_dist = pdist(pos);
pairwise_dist(pairwise_dist <= 0) = [];
pds_stats = get_stats(pairwise_dist');
vol_stats = get_stats(vol);
sur_stats = get_stats(sur);
rat_stats = get_stats(rat);
cms_stats = get_stats(cmass_dist);

stack = [pds_stats; cms_stats; vol_stats; sur_stats; rat_stats];

% Create table
total     = table(stack(:, 1), stack(:, 2), stack(:, 3), stack(:, 4),  ...
     stack(:, 5), stack(:, 6), stack(:, 7), stack(:, 8), ...
    'VariableNames', {'Mean', 'SD', 'Q05', 'Q25', 'Q50', 'Q75', 'Q95', 'Count'}, ...
    'RowNames', {'Pairwise distance', 'Mass center distance', 'Volume', ...
    'Surface area', 'Aspect ratio'});
writetable(total, path, 'WriteRowNames', true, ...
    'Sheet', [prefix, '_fine']);
end

function out = get_stats(data)
out = [...
    mean(data, 1)', ...
    std(data, [], 1)', ...
    quantile(data, [0.05, 0.25, 0.50, 0.75, 0.95], 1)', size(data, 1)];
end

