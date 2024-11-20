%
% This function can be used to screen the observer positions on <input> and their
% respective coordinates. This assumes that each column in the <input> matrix
% corresponds to a different receiver
%
% Gil Felix Greco, Braunschweig 15.09.2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function screen_receivers(input, PATH)

n_observer = size(input,2); % number of observers (columns) on input

fprintf( '\n*--------------------------------------------------------------------------*' );
fprintf( '\nLog from PANAM_SQAT_data_conversion' );
fprintf( '\n\nInput file path: %s' ,char(PATH));
fprintf( '\n\nNumber of receiver(s) found: %i \n\n',n_observer );

for i=1:n_observer
    
    fprintf( 'Receiver position %i: x = %.5g\t| y = %.5g \t| z = %.5g\n',i, input{1,i}.xobs, input{1,i}.yobs, input{1,i}.zobs );
    
end

fprintf('*--------------------------------------------------------------------------*\n' );

end