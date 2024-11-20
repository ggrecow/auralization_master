function ita_raven_createLoudspeakerSetupFile(coordinates, center)
    % Function to write source positions of the reproduction room into 
    % Speaker.ini for VBAP reproduction
    %
    % ! This version does not support lsp-directivities, as they are
    % not supported by the current VBAP version (07-Jan-2014)!
    % Input is in Raven Coordinates
    % Output is a Speaker.ini
    %
    % Authors: Sönke Pelzer, Jonas Sautter, Michael Kohnen
    % Email: mko@akustik.rwth-aachen.de
    %
    % date:     2014/06/17
    %
    % <ITA-Toolbox>
    % This file is part of the application Raven for the ITA-Toolbox. All rights reserved.
    % You can find the license for this m-file in the application folder.
    % </ITA-Toolbox>
    
    MATRIX_MATLAB2RAVEN = [0 0 -1; -1 0 0; 0 1 0];  % = MATRIX_MATLAB2SPAX
    
    radiusSpeakers=1;   %Scaling factor(!) for the distance
    %of the loudspeakers to the center point

    coordinates_mat=coordinates;
    % Convert coordinates from Raven to Matlab in two steps:
    % 1. Convert from Raven2Matlab
    %coordinates.cart = RavenProject.pSU2RVN(coordinates.cart);
    % 2. Convert from Matlab2Spax
    coordinates.cart = coordinates.cart * MATRIX_MATLAB2RAVEN;


    % set center position
    if nargin<2
        center=[0 0 0];
    end

    % open file, empty it and write
    %fid = fopen('../RavenInput/Speakers.ini', 'w');
    fid = fopen('Speakers.ini', 'w');

    % write general settings
    fprintf(fid, '\n[General]\nCenter = %d,%d,%d', center);
    fprintf(fid, '\nRadiusSpeakers = %d', radiusSpeakers);
    fprintf(fid, '\nMode = SIRERA');

    % write lsp set-up
    fprintf(fid, '\n[SIRERA]\nNumberSpeakers = ');
    array = coordinates.cart;
    fprintf(fid, '%i\nRouting = 0\nSpeakers = ', coordinates.nPoints);
    fprintf(fid, '%f , %f , %f, ', transp(array));
    tri = DelaunayTri(coordinates_mat.cart); % Calculation of the triangluar sections for VBAP
    tri = tri.convexHull();
    fprintf(fid, '\nNumberSections = %i\nSections = ', length(tri));
    fprintf(fid, '%i , %i , %i, ', transp(tri));

    %close file
    fclose(fid);
end

