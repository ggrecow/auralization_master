function Q = get_ground_reflection_coefficient( freq, sigma_e, theta, r2, soundSpeed)
% function Q = get_ground_reflection_coefficient( freq, sigma_e, theta, r2, soundSpeed)
%
%   Implementation of ground reflection coefficient considering a spherical wave correction factor.   
%   Formulas according to :
%
%   [1] Michael Arntzen (2014), Aircraft noise calculation and synthesis in a
%        non-standard atmosphere, PhD thesis, Delft university (chapter 3.2)
%
%   [2] Reto Pieren (2018), Auralization of Environmental Acoustical Sceneries - Synthesis of Road
%        Traffic, Railway and Wind Turbine Noise,  PhD thesis, Delft university (chapter 3.3)
%
%
%   This code uses the implementation of the Faddeeva functions by Johnson 
%   (because MATLAB fuctions only accepts real numbers or symbolic math):
%
%   Steven G. Johnson (2024). Faddeeva Package: complex error functions 
%   (https://www.mathworks.com/matlabcentral/fileexchange/38787-faddeeva-package-complex-error-functions), 
%    MATLAB Central File Exchange. Retrieved March 18, 2024. 
%
%
%  Description of relevant variables:
%
%   Q : ground reflection coefficient
%   Rp : plane wave reflection coefficient
%   Zn : normalized (by rho*c) ground impedance
%   F : spherical wave correction factor
%   mu : numerical distance
%   r2 : distance propagated by the ground reflection ray
%   soundSpeed : sound speed (m/s)
%
%   INPUTS:
%
%   freq : row vector [nFreq x 1]
%   frequency vector
%
%   sigma_e : scalar
%   effective flow resistance [kPa/m^2.s]
%
%   theta : scalar
%   angle between the ground plane and the ground reflected path [rad]
%
%   r2 : scalar
%   total propagation distance of the ground reflected ray [m]
%
%   OUTPUTS:
%
%   Q : row vector [nFreq x 1]
%   freq dependent ground reflection coefficient (mag & phase)
%
% Gil Felix Greco. Braunschweig 14.03.2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define ground impedance model
% model = 'Delany_Bazley';
model = 'miki';

% compute normalized ground impedance
Zn = get_ground_impedance( freq, sigma_e, model );

Zn( isnan( Zn ) ) = 0; % truncate NaN values
Zn( isinf( Zn ) ) = 0; % truncate inf values

% plane wave reflection coefficient
Rp = ( Zn.*sin(theta) - 1 ) ./ ( Zn.*sin(theta) + 1 );

% compute numerical distance
% mu = sqrt( 1j .* ( ( 2 .* pi .* freq ) ./soundSpeed ) .* ( r2 ./ 2 ) .* ( ( sin( theta ) + ( 1./Zn ) ).^2 ) ./ (1 + ( sin( theta ) ./ Zn ) ) ); % according to Arntzen, eq (3.19) 
mu = (1+1j) ./2 .* sqrt( ( 2 .* pi .* freq .* r2 ) ./soundSpeed ) .* ( sin( theta ) + 1./Zn );   % according to Pieren, eq (3.31)

mu( isnan( mu ) ) = 0; % truncate NaN values
mu( isinf( mu ) ) = 0; % truncate inf values

% compute spherical wave correction factor
% F = 1 + 1i .* mu .* sqrt( pi ) .* exp( -mu.^2 ) .* Faddeeva_erfc( -1i .* mu ); % according to Arntzen, eq (3.18) - gives a lot of inf and NaN values
F = 1 + 1i .* mu .* sqrt( pi ) .* Faddeeva_w( mu ); % according to Pieren, eq (3.30)

F( isnan( F ) ) = 0; % truncate NaN values
F( isinf( F ) ) = 0; % truncate inf values

% F = 0; % this retrieves back the plane ground reflection coefficient

% compute ground reflection coefficient
Q = Rp + ( 1 - Rp ).*F;

%% function  - compute normalized characteristic ground impedance

    function Zn = get_ground_impedance(freq, sigma_e , model)
        % function Zn = get_ground_impedance(freq, sigma_e , model)
        %
        % compute normalized (by rho*c) ground impedance model according to:
        %
        % [3] M. E. Delany and E. N. Bazley, Acoustical properties of fibrous absorbent materials,
        %      Applied Acoustics (3), 1970, pp. 105-116
        %
        % and
        %
        % [4] Y. Miki (1990), Acoustical properties of porous materials -
        %      Modifications of Delany-Bazley models, J. Acoust. Soc. Jpn
        %
        % Coefficients of those two models are quite contradictory in literature. Therefore, some values 
        % (which I believe are correct!!!) were taken from:
        %
        % [5] T. Komatsu, Improvement of the Delany-Bazley and Miki models
        %        for fibrous sound-absorbing materials, Acoust. Sci. & Tech., 29, 2 , (2008)
        %
        %
        % INPUTS:
        %
        %   freq : row vector [nFreq x 1]
        %   frequency vector
        %
        %   sigma_e : scalar
        %   effective flow resistance [kPa/m^2.s]
        %
        %   model : string
        %   apply semi empirical coefficients, according to 
        %   'Delany_Bazley' - from  Ref. [5] Table 1
        %   'miki' - from Ref. [4], Eq, (30) and (31) and Ref. [5] Table 2
        %
        %   OUTPUT:
        %
        %   Zn : row vector [nFreq x 1]
        %   normalized ground impedance as a function of freq for a given <sigma_e>
        %
        % Gil Felix Greco. Braunschweig 18.03.2024
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        switch model
            case 'Delany_Bazley'

                a = 0.0497; % Ref. [5] Table 1
                b = -0.754;
                c = 0.0758;
                d = -0.732;

                % a = 9.08; % need to multiply 1e3 to X, check https://apmr.matelys.com/PropagationModels/MotionlessSkeleton/DelanyBazleyModel.html
                % b = - 0.75;
                % c = 11.9;
                % d = - 0.73;

            case 'miki'

                a = 0.07; % from Ref. [4], Eq, (30) and (31) and Ref. [5] Table 2
                b = -0.632;
                c = 0.107;
                d = -0.632;

                % a = 5.50;  % need to multiply 1e3 to X, check https://apmr.matelys.com/PropagationModels/MotionlessSkeleton/DelanyBazleyMikiModel.html
                % b = - 0.632;
                % c = 8.43;
                % d = - 0.632;

        end

        % X = 1e3.*(freq ./ sigma_e);
        X = (freq ./ sigma_e);
        Zn =  1 + ( a .* X.^b ) - ( 1i .* c .* X.^d );

        % figure;
        % plot( freq, real(Zn) ); hold all;
        % plot( freq, imag(Zn) );
        % xlabel( 'Frequency (Hz)' );
        % ylabel( 'Normalized specific impedance');
        % legend( 'Real', 'Imaginary' );
        % set( gcf, 'color', 'w' );

    end

end