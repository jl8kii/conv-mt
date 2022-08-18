function main(arg)
if ~isdeployed
    arg = 'problem.ini';
end
opt = ini2struct(arg);
fprintf('Parsing problem definition...\n');

files_path = fullfile(opt.input, '\*.tif');
files = dir(files_path );
if numel(files) < 1
    fprintf('Failed to find files in {%s}, exiting...\n', files_path);
    return;
else
    channel = str2double(opt.parameters.channel);
    fprintf('Found %d files, reading channel %d...\n', numel(files), channel);
end
data    = {};
for i=1:numel(files)
    [~, filename, ~] = fileparts(files(i).name); 
    parts            = split(filename, '_');
    if numel(parts) <= 1, continue; end
    
    z_pos = str2double(parts{3}(2:end));
    chann = str2double(parts{4}(3:end));
    
    if chann == channel
        data{1+end} = permute(...
            im2single(imread(...
                fullfile(files(i).folder, files(i).name))), [3, 1, 2]);
    end
end
data = permute(cell2mat(data'), [2, 3, 1]);
fprintf('Sucessfully read %d images.\n', size(data, 3));

% Get working volume
xy_size = str2double(opt.parameters.xy_size);
z_size = str2double(opt.parameters.z_size);

xy_um = size(data, 1) * xy_size;
z_um  = size(data, 3) * z_size;
fprintf('Size of working area is %.1f x %.1f x %.1f μM.\n', ...
    xy_um, xy_um, z_um);

fprintf('Filtering data...\n');
data = imgaussfilt3(data, 0.25);

fprintf('Suppressing edges...\n');
sup = zeros(size(data, [1, 2]));
sup(floor(size(sup, 1)/2), floor(size(sup, 2)/2)) = 1;
sup    = bwdist(sup);
sup    = repmat(sup < median(sup(:)), [1, 1, size(data, 3)]);
data   = data .* sup;

% Calculate characteristic size of mitochondria - in that box we will
% search candidates
max_ln = str2double(opt.adjustable.max_mito_len);
fat_fc = 1.25;
xy_px  = round_odd(max_ln * fat_fc / xy_size);
z_px   = round_odd(min(max_ln * fat_fc / z_size, size(data, 3)));

% Get a reduced mesh
active = data;
x_side = single(1:size(active, 1));
y_side = single(1:size(active, 2));
z_side = single(1:size(active, 3));
[xx, yy, zz] = meshgrid(x_side, y_side, z_side);


% Growth horizontally (!) 
trials   = str2double(opt.adjustable.trial_steps);
fprintf('Iterative detection of mitochondria up to %d objects...\n', ...
    trials);
detected = [];
int      = imgaussfilt3(active, 3.0, 'FilterSize', [xy_px, xy_px, z_px]);
for steps=1:trials
    % High sum corresponds to possible mitochondria    
    [~, idx]  = max(int(:));
    [r, c, d] = ind2sub(size(active), idx);
    
    % Get a bounding box
    start     = ceil([...
        r-xy_px/2, ...
        c-xy_px/2, ...
        max(1, d-z_px/2)]);
    final     = [...
        start(1)+xy_px, ...
        start(2)+xy_px, ...
        min(start(3) + xy_px, size(active, 3))];
    
    % Do the fit here
    sub = active(...
        start(1):final(1), ...
        start(2):final(2), ...
        start(3):final(3));
    
    fit    = fit_mito(sub, xy_size, z_size);
    if isempty(fit)
        fprintf('Premature exit on %d iteration: too small object !\n', ...
            steps);
        break;
    end
    record = struct('fit', fit, 'pos', ...
        [r, c, d] .* [xy_size, xy_size, z_size], 'start', start, ...
        'final', final);
    detected = [detected; record];
    
    % Get a 3 times higher start and final
    start_fat = ceil([...
        r-3*xy_px/2, ...
        c-3*xy_px/2, ...
        max(1, d-3*z_px/2)]);
    final_fat = [...
        start_fat(1)+3*xy_px, ...
        start_fat(2)+3*xy_px, ...
        min(start_fat(3) + 3*xy_px, size(active, 3))];
    
    % Get the xyz
    xq = reshape(xx(...
        start_fat(1):final_fat(1), ...
        start_fat(2):final_fat(2), ...
        start_fat(3):final_fat(3)), [], 1) .* xy_size;
    
    yq = reshape(yy(...
        start_fat(1):final_fat(1), ...
        start_fat(2):final_fat(2), ...
        start_fat(3):final_fat(3)), [], 1) .* xy_size;
    
    zq = reshape(zz(...
        start_fat(1):final_fat(1), ...
        start_fat(2):final_fat(2), ...
        start_fat(3):final_fat(3)), [], 1) .* z_size;
    
    dq = (final_fat - start_fat) + 1;
    
    % Evaluate pdf on shifted coord    
    scaled_st = start .* [xy_size, xy_size, z_size];
    val = fit.pdf([yq - scaled_st(1), xq - scaled_st(2), zq - scaled_st(3)]);
    val = reshape(val, dq);
   
    % Normalize max to one
    val = val ./ max(val(:));
    val = 1 - val;
    
    active(...
        start_fat(1):final_fat(1), ...
        start_fat(2):final_fat(2), ...
        start_fat(3):final_fat(3)) = active( ...
        start_fat(1):final_fat(1), ...
        start_fat(2):final_fat(2), ...
        start_fat(3):final_fat(3)) .* val;   
    
     int(...
        start_fat(1):final_fat(1), ...
        start_fat(2):final_fat(2), ...
        start_fat(3):final_fat(3)) = int( ...
        start_fat(1):final_fat(1), ...
        start_fat(2):final_fat(2), ...
        start_fat(3):final_fat(3)) .* val;  
    
    if mod(steps, 50) == 0
        fprintf('Done %d from %d steps.\n', steps, trials);
    end
end
fprintf('Detecting finished, filtering...\n');
% Parse filtering criteria
min_mito_rat = str2double(opt.adjustable.min_mito_rat);
max_mito_rat = str2double(opt.adjustable.max_mito_rat);

min_mito_vol = str2double(opt.adjustable.min_mito_vol);
max_mito_vol = str2double(opt.adjustable.max_mito_vol);

min_mito_srf = str2double(opt.adjustable.min_mito_srf);
max_mito_srf = str2double(opt.adjustable.max_mito_srf);

% Add fields
for n=1:numel(detected)
    detected(n).volume = [];
    detected(n).surface = [];
    detected(n).ratio = [];
end

good = [];
for n=1:numel(detected)    
    axe = sort(eig(detected(n).fit.Sigma), 'Descend');
    
    volume  = 4*pi*prod(axe)/3;
    surface = elli_surf(axe(1), axe(2), axe(3));
    ratio   = axe(1) / min(axe);
    
    detected(n).volume  = volume;
    detected(n).surface = surface;
    detected(n).ratio   = ratio;
    
    % Check criterias
    c_ratio = ratio > min_mito_rat && ratio < max_mito_rat;
    c_vol   = volume > min_mito_vol && volume < max_mito_vol;
    c_srf   = surface > min_mito_srf && surface < max_mito_srf;
    
    if c_ratio && c_vol && c_srf
        good = [good; detected(n)];
    end
end
fprintf('%d detections deemed to be worthy, saving...\n', numel(good));
dump_result(good, opt.output, opt.name);
if ~strcmp(opt.debug_image, '0')
    image = render_debug(good, size(active), xy_size, z_size, xy_px, z_px);
    imwrite(image, opt.debug_image);
end
end

function result = elli_surf(a, b, c)
%ELLI_SURF Ellipsoid surface (approximate)
p = 1.6075;
a = a^p;
b = b^p;
c = c^p;
result = 4 * pi * ((a*b + a*c + b*c) / 3)^(1/p);
end

function S = round_odd(S)
% Round to nearest odd integer.
idx = mod(S,2)<1;
S = floor(S);
S(idx) = S(idx)+1;
end