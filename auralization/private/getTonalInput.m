
function output = getTonalInput(input, time_PANAM, source,  tag_auralization )

global pref

global input_file

global save_mat_fig

%% get tonal content over time as from PANAM

switch source
    
    case 'buzzsaw'
        %
        % Synthesis of buzzsaw tones. Procedure based on the work of:
        %   
        % Rizzi, Aumann, Lopes and Burley Auralization of Hybrid Wing–Body Aircraft Flyover
        % Noise from System Noise Predictions, JOURNAL OF AIRCRAFT Vol. 51, No. 6, 
        % DOI: 10.2514/1.C032572
        %
        %   1) fundamental freq and its <nHarmonics> are calculated based
        %       on the settings of the engine (i.e. N1 and number of
        %       blades). Because this auralization framework works with
        %       doppler-shifted EMISSION values, the doppler-effect is
        %       considered in this process.
        %
        %   2) The buzzsaw emission spectra provided by PANAM is used to obtain the SPL of 
        %        the tones calculated in step 1. The steps used for that are:
        %
        %       2.1) The 1/3-OB in which the freq of the tones calculated in step 1 lays inside 
        %               is obtained 
        %
        %       2.2) The 1/3-OB in which the freq of the emission spectra provided by PANAM    
        %               lays inside is obtained   
        %
        %       2.3) The SPL predicted within a 1/3-OB of Step 2.2 is
        %              distributed to the same the tones calculated in Step
        %              2.1. If more then one tone calculated within Step 1
        %              lays within the same 1/3-OB, the amplitude predicted
        %              by PANAM for that 1/3-OB is distributed equally to all tones  
        %
        % Gil Felix Greco, Braunschweig 18.04.2024
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % number of rotor blades
        % nBlades = str2double( input_file.n_blades );   
  
        % Max. n1 of the engine, in rotations per minute 
        maxRotationsPerMinute = str2double( input_file.max_rotations_per_minute );  

        % total number of tones to be sinthesised (i.e. 1 fundamental tone + (nHarmonics-1) harmonics)
        nHarmonics = str2double( input_file.n_harmonics );

        % Bandwidth metric, e.g., b=1 for octave filter or b=3 third octave filter.
        b = 3;  

        % get 33 octave bands from 16 Hz to 25 kHz
        [f1, ~, f2] = get_octave_bands(b);  

        % number of time steps
        nTimes = length(time_PANAM);  

        % initialize output vectors
        tonesSPLTime = zeros( nHarmonics, nTimes );  
        tonesFreqTime = zeros( nHarmonics, nTimes );
        idx_synthesised2panam = zeros( nHarmonics, nTimes );
        fundamental_freq = zeros( 1, nTimes );
        theta_doppler = zeros( 1, nTimes );
        doppler_factor = zeros( 1, nTimes );
        NfreqBSN = zeros( nHarmonics, nTimes );
        dopplerShiftedNfreqBPF = zeros( nHarmonics, nTimes );

        for i = 1:nTimes % number of source_time steps

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Step 1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % get (fundamental) Blade Passing frequency, in Hz
            fundamental_freq(i) =  ( ( input{i}.n1 ./ 100 ) * (maxRotationsPerMinute ./ 60 ) ) ; % fundamental freq of the buzzsaw tone for each time step [1 x nTimes]

            % get j-th fundamental freq and nHarmonics for each i-th time step
            NfreqBSN(1:nHarmonics,i) = (1:nHarmonics) .* fundamental_freq(i); %  fundamental freqBPF + nHarmonics for each time step [nHarmonics x  nTimes ]

            % get j-th fundamental freq and nHarmonics (doppler shifted) for each i-th time step
            % phi_x is the polar angle; phi_y is the azimuthal angle

            % theta_doppler(i) = deg2rad ( sqrt( (input{i}.phix).^2 + (input{i}.phiy).^2 ) ) ;
            theta_doppler(i) = deg2rad ( input{i}.phix )  ;

            doppler_factor(i) = 1 - ( input{i}.vel ./ input{i}.sound_speed ) .* cos ( theta_doppler(i) ) ;

            % %% law cosine approach with vectors (wrong!!! )%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % s(i,:) = [input{i}.xpos input{i}.ypos input{i}.alt];
            % 
            % r(i,:) = [input{i}.xobs input{i}.yobs input{i}.zobs];
            % 
            % doppler_factor(i) = 1 - ( input{i}.vel ./ input{i}.sound_speed ) .* cos ( deg2rad (theta_doppler(i)) ) ;

            % %% spherical to cartesian approach (wrong!!! )%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % theta_doppler(i) =  ( sin( deg2rad( input{i}.phix ) ) .* cos( deg2rad( input{i}.phiy) ) ) ;
            % 
            % doppler_factor(i) = 1 - ( input{i}.vel ./ input{i}.sound_speed ) .*  theta_doppler(i)  ;

            dopplerShiftedNfreqBPF( (1:nHarmonics) ,i ) =  NfreqBSN(1:nHarmonics,i) ./ doppler_factor(i) ;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Step 2.1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % get index of the 1/3-OB  where each tone is inside
            idx = dopplerShiftedNfreqBPF(:,i) >= f1 & dopplerShiftedNfreqBPF(:,i) <= f2; %  [nHarmonics x nBands]

            % number of tones inside each 1/3-OB
            nTonesPerToc = sum( idx~=0 );

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Step 2.2
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            a = input{i}.engine_buzzsaw(:,1)~=0; % get only tones with freq content (for receiver data, some tones may have zero SPL)
            tones(:,1) = input{i}.engine_buzzsaw(a,1); % freq vector (doppler shifted) from PANAM
            tones(:,2) = input{i}.engine_buzzsaw(a,2); % SPL (doppler shifted) from PANAM

            % get index of the 1/3-OB where each predicted buzzsaw noise is
            idx_PANAM = tones(:,1) >= f1 & tones(:,1) <= f2; % [nHarmonics x nBands]

            % number of buzzsaw noise components inside each 1/3-OB
            nTonesPerTocPANAM = sum( idx_PANAM~=0 );

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Step 2.3
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            for B = 1:nHarmonics

                if any(idx(B,:) ~= 0) % guarantee that idx(B,:) has non-zero values

                    % index where the synthesised tone lies
                    idx_synthesised =  idx(B,:) ~=0 ;  

                    % index where the synthesised tone lies in the PANAM freq vector
                    count = find( idx_PANAM( :, idx_synthesised)~=0 );

                    % when matching the amplitude of the 1/3-OB predicted
                    % by PANAM with the calculated tones, three cases may happen.

                    if isempty( count ) % 1) PANAM has no values for the synthesised tone

                        % truncate to value of this tone in the previous time-step
                        idx_synthesised2panam( B, i ) = idx_synthesised2panam( B, i-1 ); % WARNING: this may not work always, a better solution needs to be figured out

                        % get SPL (predicted from PANAM) of the 1/3-OB containing the calculated freq of the tone
                        SPL_buzzsaw = tones( idx_synthesised2panam( B, i ), 2 );

                    elseif size( count, 1) == 1 % 2) PANAM predicts amplitude of the synthesised tone in a single band

                        idx_synthesised2panam( B, i ) = count;

                        % get SPL (predicted from PANAM) of the 1/3-OB containing the calculated freq of the tone
                        SPL_buzzsaw = tones( idx_synthesised2panam( B, i ), 2 );

                    elseif size( count, 1) > 1 % 3) PANAM has 2 values predicted in the same band as the synthesised tone

                        % <idx_synthesised2panam( B, i )> would not accept a vector here, thus
                        % we fill it with the first value of <count>
                        idx_synthesised2panam( B, i ) = count(1);

                        % get energetic sum of the SPL (predicted from PANAM) in all the 1/3-OB containing the calculated freq of the tone
                        SPL_buzzsaw = 10.*log10( sum( 10.^( tones( count(1:end), 2 ) ./10 ) ) );

                    end

                    % convert from spl to pascal -  p_rms^2
                    pressure_buzzsaw = pref^2*10.^( SPL_buzzsaw /10 ); % convert from spl to pascal - % p_rms^2 vs time in current band

                    % distribute predicted (PANAM) pressure among the number of
                    % calculated tones laying within the same 1/3-OB
                    pressure_tone = sqrt( pressure_buzzsaw ./ nTonesPerToc( idx(B,:) ) );

                    % convert from p->SPL to get final amplitude over time of each <nHarmonics> tone for auralization
                    tonesSPLTime(B, i) = 20.*log10( pressure_tone./pref );

                elseif any(idx(B,:) == 0) % means idx(B,:) has only zeros (i.e. the calculated BPF tone is not inside any 1/3-OB - can happen for very low/high freq tones )

                    tonesSPLTime(B, i) = 0; % truncate to zero SPL

                end

            end

        end

        tonesFreqTime = dopplerShiftedNfreqBPF; % freq over time of each <nHarmonics> tone for auralization

        % in case of approach procedure, the buzzsaw tone is very small and
        % will only introduce more computational time to the auralization
        % procedure. However, the way that this auralization procedure is conceptuallized does 
        % not allow to have any pior knowledge which flight procedure its being auralized. 
        % Here, we try to avoid that by looking into the 1st BPF's magnitude over time, and looking  
        % if an arbitrary threshold is exceeded. If not, only one the first
        % tone is auralized, saving time during the synthesis procedure
        % (next step, using the <tonalSynthesis.m> function)

        threshold_buzzsaw = 20; % (dB SPL)
        if max( tonesSPLTime( 1, : ) ) < threshold_buzzsaw

            tonesSPLTime = tonesSPLTime ( 1, : ) ;
            tonesFreqTime = tonesFreqTime  ( 1, : ) ;

        else 
        end
        
        %% plot angle theta_doppler %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        figure('name','PROCESSING - Buzzsaw noise synthesis - doppler effect' );
        h  = gcf;
        set(h,'Units','Inches');
        pos = get(h,'Position');
        set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

        yyaxis left;
        a = plot( rad2deg (theta_doppler) ); ylabel('$\theta$ (deg)', 'Interpreter' , 'latex'); xlabel('Time, $t$ (bins)', 'Interpreter' , 'latex'); set(gcf,'color','w');
        yline(90);  ylim([0 180]);

        % plot doppler-shifted BPF
        yyaxis right; b = plot( dopplerShiftedNfreqBPF(1,:) ./ NfreqBSN(1, :) ); ylabel('Doppler-shift, $f^\prime/f_0$ (-)', 'Interpreter' , 'latex'); xlabel('Time, $t$ (bins)', 'Interpreter' , 'latex'); set(gcf,'color','w');
        ylim([0.5 1.5]); legend([a , b], '$\theta$', '$f^\prime/f_0$', 'Interpreter', 'Latex', 'Location', 'SE' );

        if isempty(tag_auralization) % if tag_auralization is empty, dont save anything
        else
            filename = strcat(tag_auralization, '_Doppler_theta_', source);
            save_pdf = 1; save_png = 0;
            export_figures( filename, save_mat_fig, save_png, save_pdf );
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% plot predicted Buzzsaw SPL (1/3-OB) and distributed SPL over calculated BPF+nHarmonics
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
        % time bin for plots - overhead position 
        [~,t_bin] = min( abs(rad2deg (theta_doppler) - 90) );  

        PLOT_buzzsaw_PANAM_synthesis(input{t_bin}.engine_buzzsaw(:,1), input{t_bin}.engine_buzzsaw(:,2), tonesFreqTime(:, t_bin), tonesSPLTime(:, t_bin), f1, f2, source, tag_auralization)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % original (deprecated approach in 17.04.2024 - takes input data from PANAM and assumes its all tones)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % j = 1;
        % for i = 1:length(time_PANAM) % number of source_time steps
        %
        %     a = input{i}.engine_buzzsaw(:,1)~=0; % get only tones with freq content (for receiver data, some tones may have zero SPL)
        %     tones(:,1) = input{i}.engine_buzzsaw(a,1);
        %     tones(:,2) = input{i}.engine_buzzsaw(a,2);
        %
        %     for nTones = 1:length(tones) % loop over the number of tones i-th source_time steps
        %         tonesFreqTime(nTones,j) = tones(nTones,1); % freq over time of each <nTones> for auralization
        %         tonesSPLTime(nTones,j) = tones(nTones,2);  % SPL over time of each <nTones> for auralization
        %     end
        %
        %     j = j+1;
        %     clear a
        % end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'fan_harmonics'

        j = 1;
        for i = 1:length(time_PANAM) % number of source_time steps

            a = input{i}.fan_harmonics(:,1)~=0; % get only tones with freq content (for receiver data, some tones may have zero SPL)
            tones(:,1) = input{i}.fan_harmonics(a,1);
            tones(:,2) = input{i}.fan_harmonics(a,2);

            for nTones = 1:length(tones) % loop over the number of tones i-th source_time steps
                tonesFreqTime(nTones,j) = tones(nTones,1); % freq over time of each <nTones> for auralization
                tonesSPLTime(nTones,j) = tones(nTones,2);  % SPL over time of each <nTones> for auralization
            end

            j = j+1;
            clear a
        end

end

%% output

output.tones = size( tonesFreqTime,1 ); % number of tones (assumes they dont change over time)
output.tonesFreqTime = tonesFreqTime; % freq of the tones over time
output.tonesSPLTime =  tonesSPLTime;  % SPl of the tones overtime

end