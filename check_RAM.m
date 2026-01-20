
% Get the information about variables in the workspace
vars_info = whos;

% Sum the Bytes used by all variables
total_bytes = sum([vars_info.bytes]);

% Convert bytes to gigabytes
total_gb = total_bytes / 1e9;

% Display the total memory in gigabytes
disp(['Total memory used by variables: ', num2str(total_gb), ' GB']);