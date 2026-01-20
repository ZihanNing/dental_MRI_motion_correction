function FOV_range = compute_FOV_range(rec,file)
% compute_FOV_range
% Read location information and compute FOV range
%
% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 17-Oct-2025

    % Define location file path
    loc_file = fullfile(rec.Names.pathOu, file);
    
    % Check if file exists
    if ~isfile(loc_file)
        error('Location file not found: %s', loc_file);
    end

    % Read numeric values from file
    fid = fopen(loc_file, 'r');
    if fid == -1
        error('Cannot open location file: %s', loc_file);
    end
    
    loc = fscanf(fid, '%f');  % read all numeric values
    fclose(fid);

    % Validate the read values
    if numel(loc) < 3
        error('Invalid location file: expected 3 numbers, got %d.', numel(loc));
    end

    loc = loc(1:3);  % ensure only first three values are used
    
    if ~(loc(1) > loc(2) && loc(2) > loc(3))
        error('Invalid location values: should satisfy num_1 > num_2 > num_3.');
    end

    % Compute FOV_range
    if ~isfield(rec, 'y') || isempty(rec.y)
        error('rec.y is missing or empty.');
    end

    FOV_range = loc ./ size(rec.y, 1);
end
