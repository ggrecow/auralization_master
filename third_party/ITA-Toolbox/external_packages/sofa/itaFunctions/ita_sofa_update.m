function ita_sofa_reinstall_current_release()
%
%   Update SOFA installation
%       delete old -> grab current
%
%   Only updates SOFA if installed as part of ita-toolbox to avoid problems
%   when messing with manual installs
%


%% check if sofa is installed
if ~exist('SOFAstart.m','file')
    ita_verbose_info('No SOFA install found, installing from scratch')
    ita_install_sofa();
    return;
end


%% check if sofa is in toolbox path
path = which('SOFAstart.m');

if ~contains(path,ita_toolbox_path)
    error(['The SOFA API seems to be installed seperately from the the ita-toolbox at %s\n',...
        'Please manually update your installation or delete it and run ita_sofa_install afterwards'],path)
end

%% remove sofa installation
ita_verbose_info('Removing old install...\n')

fullpath = fileparts(which('ita_sofa_install.m')); %taken from sofa install
[itaSofaPath,~] = fileparts(fullpath);
sofaInstallPath = fullfile(itaSofaPath,'sofa');

rmpath(genpath(sofaInstallPath))
rmdir(sofaInstallPath,'s')

ita_verbose_info('Getting a fresh version ...');
ita_sofa_install();

end %function
