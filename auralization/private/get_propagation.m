function TF = get_propagation(input, receiver, nfft, time, emission_angle_panam, show, tag_auralization)
% function TF = get_propagation(input, receiver, nfft, time, emission_angle_panam, show, tag_auralization)
%
% This function computes the transfer function between the each position of the aircraft during 
% its flight trajectory (as is from PANAM input, without any interpolation) and the receiver position.
% This transfer function corresponds to the sound attenuation imposed by the sound propagation 
% through the atmosphere. The ART ray tracing tool [1], v2023a, which can be downloaded from [2],
% is used for this purpose. It can consider a stratified atmosphere, i.e. temperature and humidity 
% vertical profiles, as well as wind, which will have in impact on the local sound speed, causing refraction. 
% The atmosphere needs to be defined here. Auralization parameters are derived from the eigenrays,
% such as geometrical spreading (distance dependent), air absorption (freq and dist dependent), 
% ground reflection (freq and distance dependent) and phase difference between direct and 1st order 
% reflected eigenrays. In the end, the transfer function of 1) direct sound path, 2) 1st order 
% reflected sound path, and 3) combination of the two aforementioned eigenrays, are
% computed. 
%
%   For more infos about the input parameters considered in the atmosphere:
%   type <help StratifiedAtmosphere>
%
%   For more infos about the propagation model used to derive the atmospheric transfer function:
%   type <help AtmosphericPropagation>
%
%   [1] Philipp Schäfer & Michael Vorländer (2021),  Atmospheric Ray Tracing: An efficient, 
%   open-source framework for finding eigenrays in a stratified, moving
%   medium, Acta Acust., 5, 26. https://doi.org/10.1051/aacus/2021018  
%
%   [2] https://www.virtualacoustics.org/GA/art/#download (Last viewed 12 March 2024)
%
%
% INPUT:
%       input : struct
%       flight profile (already trimmed) from PANAM. Contain the following fields:
%       input.x - position in x-axis, in meters [N x 1]
%       input.y - position in y-axis, in meters [N x 1]
%       input.z - altitude in meters [N x 1]
%       input.TAS -  true airspeed in meters/sec [N x 1]
%       input.thrust - thrust in kN [N x 1]

%       receiver : vector
%       [x y z] coordinates of the receiver (not moving)
%
%       nfft : scalar
%       nfft or block length used to define the freq vector used to compute
%       the atmospheric transfer function. Dictates freq discretization, and 
%       thus also the min freq that can be considered in the freq vector
%
%       time : vector
%       auralization time vector (used only for plotting)
%
%       emission_angle_panam : struct
%       contain the emission angles from PANAM. Only used to plot a
%       comparison with the emission angles from ART
%
%
%   show : logical (boolean)
%   optional parameter for figures (results) display
%   'false' (disable, default value) or 'true' (enable)
%
%   tag_auralization : string
%   containg the main path (where files should be saved) and their corresponding basic info 
%   (i.e. emission/immission & aircraft & procedure & Receiver #) to be
%   included in the name of the saved files related to auralization processes
%
% OUTPUT:
%       TF : struct
%       ART/ITA results containing the atmosperic transfer function for each combination of
%       source x receiver positions (row vector). For each results, the 1st, 2nd and 3rd columns 
%       corresponds to the transfer function of the direct path, 1st order reflection, and the 
%       combination of both eigenrays, respectively.  
%
% Author: Gil Felix Greco, Braunschweig 08.12.2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global input_file
global fs % sampling frequency used for auralization, defined in the <auralization_master>
global save_mat_fig

%% Define source and receiver position

% interpolate trajectory to get more points (not clear to me if this approach is valid because we need to interpolate data in more than 1 dimension)
% interpFactor = 5;
% interpNum = size(input.x,1)*interpFactor;
% xx = linspace( input.x(1), input.x(end), interpNum);
% xx = transpose(xx);
% 
% zz = interp1( input.x, input.z, xx , 'spline' );
% yy = zeros(interpNum,1);
% 
% source = [xx yy zz];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 15.11.2024 - truncate receiver to a finite heigth
% if receiver is flush mounted in the ground, ray-tracing will have trouble 
% to find reflections. Gil have tested many options and increasing the ray
% tracing parameters (i.e. art.maxReceiverRadius and integrationTimeStep) 
% does not really helps while substantially increases the computational time.
% Truncating the receiver height to 1cm already overcomes the problem withouth
% having to change the code too much. Another alternative would be to
% double the direct path, but then what shall be done with the reflection
% coefficient?
if receiver(1,3)==0  % if receiver heigth is 0 m
    receiver(1,3) = 0.01; % truncate to 1 cm
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get trajectory from input data as is
source = [input.x input.y input.z];

% truncate trajectory to 1 point just for testing
% source = [input.x(1) input.y(1) input.z(1)];
% receiver = [2000 0 20];

%% Define atmosphere

atmos = StratifiedAtmosphere; % declare StratifiedAtmosphere class

% wind settings
atmos.windProfile = input_file.wind_profile;          % String: 'zero', 'constant', 'log'

% atmos.constWindVelocity = str2double( input_file.const_wind_velocity );         % Wind velocity for Const Wind Profile [m/s]
% atmos.constWindDirection = str2double( input_file.const_wind_direction );  % Normal in wind direction [x y z]
% atmos.surfaceRoughness = str2double( input_file.surface_roughness );        % Surface Roughness for Log Wind Profile [m]
% atmos.frictionVelocity = str2double( input_file.friction_velocity );             % Friction velocity for Log Wind Profile [m/s]

% temperature settings
atmos.temperatureProfile =  input_file.temperature_profile;   % String: 'constant', 'isa'
temperatureCelsius = str2double( input_file.temperature_celsius ); % Temperature in degree Celsius
atmos.constTemperature = temperatureCelsius + 273.15;    % Temperature for Const mode [K]

atmos.constStaticPressure = str2double( input_file.const_static_pressure ); % Static pressure used in a constant temperature profile [Pa]

% humidity settings
atmos.humidityProfile = input_file.humidity_profile;  % String: 'constant'
atmos.constRelHumidity = str2double( input_file.const_rel_humidity );        % Constant Realitive Humidity [%]

%% Define propagation model

% Initializing Propagation Model
propagationModel = AtmosphericPropagation;

% Atmosphere - It is important to use the same atmosphere as for the eigenray simulation.
propagationModel.atmosphere = atmos;

% define freq vector - add 1 bin to freq vector
nBins = ceil( (nfft+1)/2) ;     % calculate the number of unique fft points
freq = linspace(0, fs/2, nBins); % freq vector
propagationModel.frequencyVector = freq;

%% ground reflection
% The ground reflection factor is required for the transfer function
% calculation. It is 1 by default. It can either be a singe value or a
% complex-valued vector with same length as the frequency vector.

% propagationModel.groundReflectionFactor = 0.9 * ones(size(propagationModel.frequencyVector));

% sigma_e_dict = {
%     'snow': 25.,
%     'forest': 50.,
%     'grass': 250.,
%     'dirt_roadside': 500.,
%     'dirt': 5000,
%     'asphalt': 10000,
%     'concrete': 50000
%     }

% effective flow resistance [kPa/m^2.s]
sigma_e = str2double( input_file.sigma_e );  

% coarse approximation of sound speed as a function of temperature for ground reflection calculation
soundSpeed = 331.3 + (temperatureCelsius * 0.606);  % This still needs to be tested for the case when temperature profiles are considered.

%% Calculate Eigenrays and tranfer function
tic;

art = AtmosphericRayTracer; % declare AtmosphericRayTracer class

art.maxReceiverRadius = 0.1; % Maximum value for receiver radius [m]
% art.integrationTimeStep = 0.01;
% art.maxAngleForGeomSpreading = 0.001;    %Maximimum delta angle of initial direction of neighboring rays used for the calculation of the spreading loss [°]
        
TF = cell (size(source,1),1); % pre allocate memory
thetaReflectedRay = zeros(size(source,1),1);
propDistanceReflectedRay = zeros(size(source,1),1);
launchAngle_direct = zeros(1, size(source,1));
launchAngle_reflected = zeros(1, size(source,1));

for i = 1:size(source,1)
    
    % Calculate Eigenrays
    eigenrays(i,:) = art.FindEigenrays(atmos, source(i,:), receiver);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ground reflection
    % The ground reflection factor is required for the transfer function
    % calculation. It is 1 by default. It can either be a singe value or a
    % complex-valued vector with same length as the frequency vector.

    % find idx corresponding to reflection coordinate (i.e. point where reflected ray meets the ground at z = 0)
    idx_reflection = find( eigenrays(i, 2).r.z == 0 );  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (this will lead to reflected angle)

    % % distance between reflection (at z=0) and receiver, i.e. hypotenuse  
    % hypotenuse = norm( receiver - eigenrays(i, 2).r.cart(idx_reflection,:) );  
    % 
    %  % angle theta between the ground plane and the ground reflected path [rad]
    % thetaReflectedRay(i) = asin( norm( receiver ) / hypotenuse ); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (this will lead to incident angle)
    
    % distance between reflection (at z=0) and receiver, i.e. hypotenuse 
    % hypotenuse = norm( source(i,:) - eigenrays(i, 2).r.cart(idx_reflection,:) );  
    hypotenuse = norm( eigenrays(i, 2).r.cart(idx_reflection,:) - source(i,:)  );  

    % angle theta between the ground plane and the ground reflected path [rad]
    % thetaReflectedRay(i) = asin( norm( source(i,:) ) / hypotenuse ); 
    thetaReflectedRay(i) = asin( norm( source(i,3) ) / hypotenuse ); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % total propagation distance of the ground reflected ray [m]
    propDistanceReflectedRay(i) =  eigenrays(i, 2).pathLength(); 

    % launch angle - original data convention (spherical coordinate) is: 0° points upwards (northpole) and 180° downwards (southpole). Here we use: 0° (east) and 180° (west)
%     launchAngle_direct(i) = eigenrays(i,1).n0.theta_deg - 90;
%     launchAngle_reflected(i) = eigenrays(i,2).n0.theta_deg -90 ;

    % launch angle - original data convention (polar coordinate) is (i think): 0° (west) and 180° downwards (east), increasing counterclockwise. Here we use: 0° (east) and 180° (west), increasing clockwise
    launchAngle_direct(i) = 360 - eigenrays(i,1).n0.alpha_deg;
    launchAngle_reflected(i) = 360 - eigenrays(i,2).n0.alpha_deg;


    propagationModel.groundReflectionFactor = get_ground_reflection_coefficient( propagationModel.frequencyVector,... % freq (row) vector
                                                                                               sigma_e, ... % effective flow resistance [kPa/m^2.s]
                                                                                               abs( thetaReflectedRay(i) ) ,... % angle theta between the ground plane and the ground reflected path [rad]
                                                                                               propDistanceReflectedRay(i), ...   % total propagation distance of the ground reflected ray [m]
                                                                                               soundSpeed );

    % Transfer function
    % The transfer function (TF) is calculated for each eigenray individually and
    % then combined using the priciple of superposition.
    % It considers spreading loss, propagation delay, air attenuation and a
    % reflection factor for the ground reflection.

    % Calculate transfer function
    [ combinedTF, individualTFs ] = propagationModel.TransferFunction( eigenrays(i,:) );
    TF{i,1} =  ita_merge( individualTFs, combinedTF );
    
end

clear idx_reflection hypotenuse;

fprintf('\n- Propagation (ART) - Calculation time (eigenrays) was\t%f sec\n', toc);

%% plots

if show == 1
        
    % time vector
    xx = time;

    % eigenrays.plot(); % ART function - blackbox plot function from ART
    % PLOT_eigenrays(receiver, source, eigenrays); % self-programmed plot eigenrays function
    
    %% plot TF "spectrogram"
    
    % spectrogram

    % att_direct = zeros( size(freq,1), size(source,1) );
    % att_reflected = zeros( size(freq,1), size(source,1) );
    % att_combined = zeros( size(freq,1), size(source,1) );

    % for i = 1:size(source,1) % convert to dB
    %     att_direct(:,i) = 10.*log10( abs( TF{i,1}.freq(:,1) ).^2.);
    %     att_reflected(:,i) = 10.*log10( abs( TF{i,1}.freq(:,2) ).^2 );
    %     att_combined(:,i) = 10.*log10( abs( TF{i,1}.freq(:,3) ).^2 );
    % 
    %     % pref = 2e-5;
    %     % att_direct(:,i) = 10.*log10( abs( TF{i,1}.freq(:,1) ).^2./pref^2 );
    %     % att_reflected(:,i) = 10.*log10( abs( TF{i,1}.freq(:,2) ).^2./pref^2 );
    %     % att_combined(:,i) = 10.*log10( abs( TF{i,1}.freq(:,3) ).^2 ./pref^2 );
    % end

    % tag_title =  'OUTPUT - transfer function (direct path)';
    % PLOT_transfer_function(xx, freq, att_direct, tag_title)
    % 
    % tag_title =  'OUTPUT - transfer function (1st order reflection)';
    % PLOT_transfer_function(xx, freq, att_reflected, tag_title)
    % 
    % tag_title =  'OUTPUT - transfer function (direct path + 1st order reflection)';
    % PLOT_transfer_function(xx, freq, att_combined, tag_title)

    %% plot theta

    h  = figure;
    set(h, 'name', 'PROCESSING - 1st order reflection from ART (theta and propagation distance)' );
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

    yyaxis left
    plot( xx, rad2deg (thetaReflectedRay) );
    ylabel( 'Incidence angle, $\theta_{\mathrm{in}}$ (deg)', 'Interpreter', 'Latex' );
    set( gcf,'color','w' );

    % plot propagated distance reflected ray 
    yyaxis right
    plot( xx, propDistanceReflectedRay );
    xlabel( 'Time, $t$ (s)','Interpreter','Latex' );
    ylabel( 'Propagated distance, $r_2$ (m)','Interpreter','Latex' );
    set( gcf,'color','w');

    if isempty(tag_auralization) % if tag_auralization is empty, dont save anything
    else
        filename = strcat(tag_auralization, '_propagation_reflected_ray_parameters');
        save_pdf = 1; save_png = 0;
        export_figures( filename, save_mat_fig, save_png, save_pdf );
    end

    %% plot theta launch

    h  = figure;
    set(h, 'name', 'PROCESSING - launch direction from ART' );
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

    yyaxis left
    a = plot( xx, emission_angle_panam.phix, '-'); hold on;
    b = plot( xx, launchAngle_direct, '--x' ); hold on;
    % plot( xx, launchAngle_reflected,'--');
    ylabel( 'Emission angle, $\alpha^{*}$ (deg)', 'Interpreter', 'Latex' );

    b.MarkerIndices = 1:10:length(launchAngle_direct);

    % diference between panam and art (direct path)
    yyaxis right
    plot( xx, launchAngle_direct - emission_angle_panam.phix );
    xlabel('Time, $t$ (s)','Interpreter','Latex' );
    ylabel('$\alpha^{*}_\mathrm{ART}- \alpha^{*}_\mathrm{PANAM}$ (deg)','Interpreter','Latex' );
    legend([a,b], {'PANAM', 'ART - direct path'}, 'Location', 'SE' );
    legend box off;

    set( gcf,'color','w');

    if isempty(tag_auralization) % if tag_auralization is empty, dont save anything
    else
        filename = strcat(tag_auralization, '_ART_emission_angle');
        save_pdf = 1; save_png = 0;
        export_figures( filename, save_mat_fig, save_png, save_pdf );
    end

    %% Plot freq vs. 1 (arbitrarly chosen) time step using ita toolbox 
    
    % overhead_idx = round(length(TF)/2); % one source/receiver combination
    % [~, overhead_idx] =  min(propDistanceReflectedRay); % one source/receiver combination (overhead aircraft position)
    [~, overhead_idx] =  max(rad2deg (thetaReflectedRay)); % one source/receiver combination (overhead aircraft position)

    % % plot using ita toolbox
    % allTFs = TF{ overhead_idx ,1 }; % ita_merge(individualTFs, combinedTF);
    % allTFs.channelNames = {'Direct path TF', 'Ground reflection TF', 'Combined TF'};
    % allTFs.pf();
    
    att_direct = 10.*log10( abs( TF{overhead_idx,1}.freq(:,1) ).^2. );
    att_reflected = 10.*log10( abs( TF{overhead_idx,1}.freq(:,2) ).^2 );
    att_combined = 10.*log10( abs( TF{overhead_idx,1}.freq(:,3) ).^2 );

    % pref = 2e-5; % reference sound pressure
    % att_direct = 10.*log10( abs( TF{overhead_idx,1}.freq(:,1) ).^2./pref^2 );
    % att_reflected = 10.*log10( abs( TF{overhead_idx,1}.freq(:,2) ).^2./pref^2 );
    % att_combined = 10.*log10( abs( TF{overhead_idx,1}.freq(:,3) ).^2 ./pref^2 );

    h  = figure;
    set(h, 'name', 'OUTPUT - atmospheric transfer function' );
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

    a = semilogx( freq, att_combined, 'Color', [0.93,0.69,0.13], 'LineWidth', 0.5 ) ; hold on;
    b = semilogx( freq, att_direct , 'Color', [0.00,0.45,0.74], 'LineWidth', 1 ) ; 
    c = semilogx( freq, att_reflected , 'Color', [0.85,0.33,0.10], 'LineWidth', 1 ) ;  

    xlabel( 'Frequency, $f$ (Hz)', 'Interpreter', 'Latex' );

    ylabel('$10\cdot\log_{10}\left(|H_{\mathrm{atm,ART}}|^2 \right)$','Interpreter','Latex'); 
    % ylabel( 'Transmission loss, TL (dB re 20$~\mu$Pa) ', 'Interpreter', 'Latex' );

    legend([ b, c , a ], 'Direct path', 'Reflected path', 'Direct path + reflected path', 'Location', 'SW');

    xlim( [ 20 25e3 ] );
    ylim( [ -200 0 ] );

    xticks([ 20 1e2 1e3 10e3 20e3]);
    xticklabels({'20','100','1k', '10k', '20k'});

    h = gca;
    h.YAxis.MinorTick = 'on'; % Must turn on minor ticks if they are off
    h.YAxis.MinorTickValues = -200:10:0; % Minor ticks which don't line up with majors

    set( gcf, 'color', 'w' );

    if isempty(tag_auralization) % if tag_auralization is empty, dont save anything
    else
        filename = strcat(tag_auralization, '_atmospheric_transfer_function_overhead');
        save_pdf = 1; save_png = 0;
        export_figures( filename, save_mat_fig, save_png, save_pdf );
    end

else
end

end