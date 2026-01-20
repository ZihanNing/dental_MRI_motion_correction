

function [] = dat2ISMRMRD(fileName, fileNameOu, pathOu, conversionParam)

[pathTemp,file] = fileparts(fileName);
fprintf('=====  Converting .dat file to a .h5 file====\n   %s\n', file);

if nargin<2 || isempty(fileNameOu);fileNameOu=fileName;end
if nargin<3 || isempty(pathOu);pathOu=[];end
if nargin<4 || isempty(conversionParam);conversionParam='-M';end

%%% NAME HANDLING
[pathOuTemp,fileOu] = fileparts(fileNameOu);

if ~isempty(pathOu); pathOu = pathOuTemp; end
if isempty(pathOu); pathOu = pathTemp; end

%%% CONVERT
cmd = sprintf( 'siemens_to_ismrmrd -f %s.dat -o %s.h5 %s', fileName, fullfile(pathOu,fileOu), conversionParam);
dos(cmd);
 

% https://github.com/ismrmrd/siemens_to_ismrmrd 
%  This is the file that is being converted. The user must supply it using the option -f. Siemens dat file can be a Numaris/4 VB/VD/VE or Numaris/X file. In case of a file that contains multiple measurements, all measurements can be exported using the -Z option (default off) into separate files appended by the measurement number. In this example, multi-measurement file meas_MID00832.dat file has 2 measurements and two output files resulting_file_1.h5 and resulting_file_2.h5 are created:
% 
% $ siemens_to_ismrmrd -f meas_MID00832.dat -o resulting_file.h5 -Z
% The "all measurements" -Z can be combined with the -M to create a multi-measurement output file where different measurements are stored in separate groups. In this example, multi-measurement file meas_MID00832.dat file has 2 measurements, stored in HDF5 groups dataset_1 and dataset_2 respectively in output file resulting_file.h5:
% 
% $ siemens_to_ismrmrd -f meas_MID00832.dat -o resulting_file.h5 -Z -M
% A single measurement can be specified using the option -z (default value 1). Primary measurement data is usually stored as the last measurement, with earlier measurements being dependencies. In this example, the second measurement from the meas_MID00832.dat file is converted and stored in the resulting_file.h5 file:
% 
% $ siemens_to_ismrmrd -f meas_MID00832.dat -o resulting_file.h5 -z 2