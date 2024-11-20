function alpha = airAttenuationBASS(atmos, altitude, f)
%Calculates the sound attenutation for this atmosphere at a given altitude
%and frequency bins according to ISO 9613-1.
%
%Inputs:
%altitude:  Altitude (m) [1x1 double]
%f:         Frequency vector (Hz) [1xN / Nx1 double]
%
%Outputs:
%alpha:     Attenuation coefficient (dB/m)

% -------------------------------------------------------------------------
%                        ____  __________  _______
%                       //  / //__   ___/ //  _   |
%                      //  /    //  /    //  /_|  |
%                     //  /    //  /    //  ___   |
%                    //__/    //__/    //__/   |__|
%
% -------------------------------------------------------------------------
%                  ARTMatlab - ITAGeometricalAcoustics
%        (c) Copyright Institute of Technical Acoustics (ITA)
%  This file is part of the ARTMatlab application. Some rights reserved.
% You can find the license in the LICENSE.md file in the ARTMatlab folder.
%--------------------------------------------------------------------------


T = atmos.T(altitude);
hr = atmos.humidity(altitude);
pa = atmos.staticPressure(altitude);

pr = 101325;
T0 = 293.15;
T01 = 273.16;

psat = pr * 10^(-6.8346 * (T01 / T)^1.261 + 4.6151);
h = hr * (psat / pa);

frO = fRelaxO(pa, pr, h);
frN = fRelaxN(pa, pr, h, T, T0);

alpha = attenuationCoeff(f, pa, pr, T, T0, frN, frO);

function frO = fRelaxO(pa, pr, h)
%pa:    Ambient atmospheric pressure (Pa)
%pr:    Reference atmospheric pressure (101325 Pa)
%h:     Molar concentration of water vapour

frO = pa/pr*( 24 + 4.04e4*h * (0.02+h)/(0.391+h) );

function frN = fRelaxN(pa, pr, h, T, T0)
%pa:    Ambient atmospheric pressure (Pa)
%pr:    Reference atmospheric pressure (101325 Pa)
%h:     Molar concentration of water vapour
%T:     Ampient temperature (K)
%T0:    Reference temperature (293.15K)

frN = pa/pr*(T0/T).^(1/2) * ( 9 + 280*h*exp( -4.17*( (T0/T).^(1/3) - 1 ) ) );

function A = attenuationCoeff(f, pa, pr, T, T0, frN, frO)
%f:     Frequency vector
%pa:    Ambient atmospheric pressure (Pa)
%pr:    Reference atmospheric pressure (101325 Pa)
%T:     Ampient temperature (K)
%T0:    Reference temperature (293.15K)
%frN:   Relaxation frequency of nitrogen
%frO:   Relaxation frequency of oxygen

term_1 = 1.84e-11 ./ ( (T0./T).^(1/2) .* (pa./ pr) );
term_2 = (.1068.*exp(-3352./ T).*frN)./(f.^2 + frN.^2);
term_3 = (.01278.*exp(-2239.1./T).*frO)./(f.^2 + frO.^2);

A = f.^2.* ( term_1 + (T0./T).^(2.5) * (term_2 + term_3) );

% 	!-------------------------------------------------------------
% 	!Umrechnung von Neper/meter in dB/m : multiplizieren mit 20.0*log10(exp(1.0))
% 	!-------------------------------------------------------------
A = A*20*log10( exp(1.0) ); % !--> Faktor 8.686
