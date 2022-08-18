function result = fit_mito(data3, xy_size, z_size)
%FIT_MITO Meow
data3     = imgaussfilt3(data3, 0.5);
t         = graythresh(data3);
bin       = bwlabeln(data3 > t);
regs      = regionprops3(bin, 'Volume');
[~, I]    = max([regs.Volume]);
idx       = find(bin == I);

% Degenerated case
if numel(idx) < 3^5
    result = [];
    return;
end

val       = data3(idx);
val       = ceil(10 * rescale(val))+1;
[X, Y, Z] = ind2sub(size(data3), idx);
data      = [X, Y, Z];
% Scale it to micrometers
data      = data .* [xy_size, xy_size, z_size];


% Replace val times the data
arr       = cell(size(data, 1), 1);
for m=1:size(data, 1)
    arr{m} = repmat(data(m, :), [val(m), 1]);
end
arr  = cell2mat(arr);
dist = fitgmdist(arr, 1);

% Well, do not forget that its scaled !
result = dist;
end

