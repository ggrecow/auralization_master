% ISO 9613-1 Table 1: Air Temperature -20 degree Celsius

p_0 = 101.325; % kPa

f = 50; % Hz
T = -20; % degree Celsius
rh = 10; % Pecentage
t1_tneg20_example1 = ita_atmospheric_absorption_iso9613( f, T, rh, p_0 ); % dB/km
fprintf( 'Calculated %.3f, tabulated value %.3f\n', t1_tneg20_example1 * 1e3, 5.89e-1 )

% ISO 9613-1 Table 1: Air Temperature 20 degree Celsius

f = 6300; % Hz
T = 20; % degree Celsius
rh = 15; % Pecentage
t1_t20_col2_example1 = ita_atmospheric_absorption_iso9613( f, T, rh, p_0 ); % dB/km
fprintf( 'Calculated %.3f, tabulated value %.3f\n', t1_t20_col2_example1 * 1e3, 175 )

f = 50; % Hz
T = 20; % degree Celsius
rh = 80; % Pecentage
t1_t20_col2_example2 = ita_atmospheric_absorption_iso9613( f, T, rh, p_0 ); % dB/km
fprintf( 'Calculated %.3f, tabulated value %.3f\n', t1_t20_col2_example2 * 1e3, 5.01e-2 )

f = 1000; % Hz
T = 20; % degree Celsius
rh = 80; % Pecentage
t1_t20_col2_example3 = ita_atmospheric_absorption_iso9613( f, T, rh, p_0 ); % dB/km
fprintf( 'Calculated %.3f, tabulated value %.3f\n', t1_t20_col2_example3 * 1e3, 5.15 )

f = 10000; % Hz
T = 20; % degree Celsius
rh = 80; % Pecentage
t1_t20_col2_example3 = ita_atmospheric_absorption_iso9613( f, T, rh, p_0 ); % dB/km
fprintf( 'Calculated %.3f, tabulated value %.3f\n', t1_t20_col2_example3 * 1e3, 1.05e2 )
