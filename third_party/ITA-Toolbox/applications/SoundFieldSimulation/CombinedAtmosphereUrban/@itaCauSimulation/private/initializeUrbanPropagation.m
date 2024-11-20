function [urbanPropagation] = initializeUrbanPropagation()
%INITIALIZEURBANPROPAGATION 

urbanPropagation = itaPigeonProject;
urbanPropagation.run_quiet = true;

% if ~exist(urbanPropagation.result_file_path,'dir')
%     mkdir(urbanPropagation.result_file_path)
% end
% urbanPropagation.outFilePath = urbanPropagation.pidgeonPath;

end

