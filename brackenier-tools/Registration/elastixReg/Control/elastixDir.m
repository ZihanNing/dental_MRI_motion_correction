
function elastixdir = elastixDir()

version = 'elastix-4.9.0-linux';%version 5.0 did not work with libraries installed on the current machine
addpath(genpath(strcat('/home/ybr19/Software/Elastix/', version)));

elastixdir = sprintf('/home/ybr19/Software/Elastix/%s/bin/', version);    % elastix (http://elastix.isi.uu.nl/) 
elastixlib = sprintf('/home/ybr19/Software/Elastix/%s/lib/', version); 

if ispc
    setenv('PATH', elastixdir);
elseif ismac
    elastixdir = sprintf('export DYLD_LIBRARY_PATH=%s; %s',elastixlib,elastixdir);
elseif isunix
%      tt1 = sprintf('export PATH=%s;',elastixdir);
%      tt2 = sprintf('export LD_LIBRARY_PATH=%s;',elastixlib);
%      elastixdir = strcat( tt1, tt2);
     elastixdir = sprintf('export LD_LIBRARY_PATH=%s; %s',elastixlib,elastixdir);
end