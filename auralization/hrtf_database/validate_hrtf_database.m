% validate_hrtf_database.m
%
% Verification of the HUTUBS HRTF set (<hrtf_database> = HUTUBS) against the shipped
% FABIAN set. Needs the HUTUBS files of the subjects below (see <get_hrtf_database>).
%
%   1) angle convention: for sources at the left (az = 90), right (270), front (0) and
%      above (el = 90), the ear that is closer to the source must receive the sound
%      earlier and with more energy (both databases), and the frontal and overhead
%      directions must be nearly symmetric
%   2) same head in both databases: HUTUBS pp1 (FABIAN measured in the HUTUBS system)
%      vs. FABIAN, absolute third-octave band differences in the horizontal plane
%      (10 deg steps), median across directions and both ears. Brinkmann et al. (2019),
%      Fig. 7, report below 1 dB up to 8 kHz and up to 3 dB at 20 kHz (auditory filters,
%      measured HRIRs; here third-octave bands and SH-interpolated HRIRs, SH order 16 vs.
%      35, so slightly larger values are expected). A swapped or mirrored angle
%      convention gives much larger differences (shown for comparison)
%   3) diffuse-field (DF) method of <get_hrtf_df_equalization> applied to FABIAN vs. the
%      CTF provided with the FABIAN database (third-octave band levels)
%   4) DF filters of the HUTUBS subjects: deviation of the FIR filter (nTaps) from its
%      target magnitude, and comparison with the FABIAN CTF
%   5) equalization filters are cached per HRTF set (different subjects give different
%      filters)
%
% The FF equalization and the conversion between equalizations are checked for any HRTF
% set with <validate_hrtf_equalization> (set <hrtf> there).
%
% Results (Octave 8.4, 19.11.2026; HUTUBS pp1, pp22, pp26), to be confirmed in MATLAB:
% 1) passed for all sets; 2) median |dL| <= 0.4 dB up to 1 kHz, <= 1.3 dB up to 8 kHz,
% 2.8 dB at 16 kHz, 5.1 dB at 20 kHz (mirrored azimuth: 16.9 dB up to 8 kHz);
% 3) max |dL| = 0.23 dB (50 Hz - 16 kHz), -0.8 dB at 20 kHz; 4) FIR vs. target: max
% |dL| = 0.01 dB; DF of pp1 vs. FABIAN CTF: max 0.4 dB up to 8 kHz; 5) passed.
%
% Gil Felix Greco, Braunschweig 19.11.2026
%
%   AI disclosure: Generative AI assistance (Anthropic Claude) was used to develop this
%   code under the author's direction. All codes were verified by the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

subjects = [1 22 26];   % HUTUBS subjects to check (pp1 is needed for step 2)

hF = get_hrtf_database( 'FABIAN' );
for k = 1:numel( subjects )
    hH(k) = get_hrtf_database( 'HUTUBS', subjects(k) ); %#ok<SAGROW>
end
fs = hF.fs;   % 44.1 kHz, both databases

N   = 2^14;
f   = (0:N/2)' * fs/N;
fc  = 1000 * 2.^( (-13:13)/3 );            % third-octave bands 50 Hz - 20 kHz
bandLevel = @(P) cell2mat( arrayfun( @(b) 10*log10( mean( P( f >= fc(b)*2^(-1/6) & f < fc(b)*2^(1/6), : ), 1 ) ), ...
                                     (1:numel(fc))', 'UniformOutput', false ) );    % [nBands x nCols]
spec = @(h) abs( fft( h, N ) ).^2;         % power spectrum, [N x nCols]
half = @(P) P( 1:N/2+1, : );

%% 1) angle convention
fprintf( '1) angle convention (closer ear earlier and louder; front/overhead symmetric)\n' );
sets = [ {hF} num2cell( hH ) ];
ok1  = true;
for s = 1:numel( sets )
    [l, r] = get_hrir( [90 270 0 0], [0 0 0 90], sets{s} );
    E   = 10*log10( [ sum( l.^2 ); sum( r.^2 ) ] );            % [2 x 4] energy, dB
    toa = [ il_onset( l ); il_onset( r ) ];                       % [2 x 4] onset, samples
    c1  = E(1,1) > E(2,1) + 5 && toa(1,1) < toa(2,1) - 10;        % az = 90: left ear
    c2  = E(2,2) > E(1,2) + 5 && toa(2,2) < toa(1,2) - 10;        % az = 270: right ear
    c3  = abs( diff( E(:,3) ) ) < 2 && abs( diff( toa(:,3) ) ) < 3; % front: symmetric
    c4  = abs( diff( E(:,4) ) ) < 2 && abs( diff( toa(:,4) ) ) < 3; % above: symmetric
    fprintf( ['   %-12s az 90: ILD %+5.1f dB, ITD %+5.1f samples | az 270: ILD %+5.1f dB, ITD %+5.1f | ' ...
              'front: ILD %+4.1f, ITD %+4.1f | above: ILD %+4.1f, ITD %+4.1f -> %s\n'], sets{s}.name, ...
             E(1,1)-E(2,1), toa(2,1)-toa(1,1), E(1,2)-E(2,2), toa(2,2)-toa(1,2), ...
             E(1,3)-E(2,3), toa(2,3)-toa(1,3), E(1,4)-E(2,4), toa(2,4)-toa(1,4), il_pass( c1 && c2 && c3 && c4 ) );
    ok1 = ok1 && c1 && c2 && c3 && c4;
end

%% 2) same head in both databases: HUTUBS pp1 vs. FABIAN, horizontal plane
idx1 = find( subjects == 1, 1 );
if isempty( idx1 )
    fprintf( '\n2) skipped (HUTUBS pp1 not in <subjects>)\n' );
else
    az = (0:10:350)';  el = zeros( size( az ) );
    [lF, rF] = get_hrir( az, el, hF );
    [lH, rH] = get_hrir( az, el, hH(idx1) );
    dL = abs( [ bandLevel( half( spec( lH ) ) ) - bandLevel( half( spec( lF ) ) ), ...
                bandLevel( half( spec( rH ) ) ) - bandLevel( half( spec( rF ) ) ) ] );   % [nBands x 2*nDir]
    d2 = median( dL, 2 );
    fprintf( '\n2) HUTUBS pp1 vs. FABIAN, horizontal plane, median |dL| across directions and ears\n' );
    fprintf( '   %8s: %s\n', 'f (Hz)', sprintf( '%7.0f', fc ) );
    fprintf( '   %8s: %s\n', '|dL| dB', sprintf( '%7.2f', d2 ) );
    % same comparison with a mirrored azimuth (az -> -az) for HUTUBS: must be much worse
    [lM, rM] = get_hrir( -az, el, hH(idx1) );
    dM = abs( [ bandLevel( half( spec( lM ) ) ) - bandLevel( half( spec( lF ) ) ), ...
                bandLevel( half( spec( rM ) ) ) - bandLevel( half( spec( rF ) ) ) ] );
    d2max = max( d2( fc <= 8000 ) );
    dMmax = max( median( dM( fc <= 8000, : ), 2 ) );
    ok2   = d2max < 1.5 && d2max < dMmax/5;
    fprintf( ['   max up to 8 kHz: %.2f dB (expected about 1 dB, Brinkmann et al. 2019); ' ...
              'mirrored azimuth: %.2f dB -> %s\n'], d2max, dMmax, il_pass( ok2 ) );
end

%% 3) DF method applied to FABIAN vs. the CTF of the FABIAN database
EQ_ctf = get_hrtf_equalization_filter( 'DF', hF );     % CTF provided with the database
DF_F   = get_hrtf_df_equalization( hF );               % method used for HUTUBS
L_ctf  = bandLevel( half( spec( EQ_ctf.left ) ) );
L_dfF  = bandLevel( half( spec( DF_F.left ) ) );
in16   = fc <= 16000;
d3     = L_dfF - L_ctf;
ok3    = max( abs( d3( in16 ) ) ) < 0.3;
fprintf( '\n3) DF method (get_hrtf_df_equalization) on FABIAN vs. FABIAN CTF (database)\n' );
fprintf( '   %8s: %s\n', 'f (Hz)', sprintf( '%7.0f', fc ) );
fprintf( '   %8s: %s\n', 'dL dB', sprintf( '%7.2f', d3 ) );
fprintf( '   max |dL| 50 Hz - 16 kHz: %.2f dB (expected < 0.3 dB) -> %s\n', max( abs( d3( in16 ) ) ), il_pass( ok3 ) );

%% 4) DF filters of the HUTUBS subjects
fprintf( '\n4) DF filters of the HUTUBS subjects\n' );
ok4 = true;
for k = 1:numel( hH )
    DF   = get_hrtf_df_equalization( hH(k) );
    Gt   = interp1( DF.f, DF.G.^2, f, 'linear' );          % target power on the analysis grid
    d4   = bandLevel( half( spec( DF.left ) ) ) - bandLevel( Gt );
    dCTF = bandLevel( half( spec( DF.left ) ) ) - L_ctf;
    ok   = max( abs( d4( in16 ) ) ) < 0.1;
    ok4  = ok4 && ok;
    fprintf( ['   %-12s FIR (%d taps) vs. target: max |dL| = %.3f dB (50 Hz - 16 kHz) -> %s; ' ...
              'vs. FABIAN CTF: max |dL| = %.1f dB up to 8 kHz, %.1f dB at 16 kHz\n'], hH(k).name, ...
             numel( DF.left ), max( abs( d4( in16 ) ) ), il_pass( ok ), ...
             max( abs( dCTF( fc <= 8000 ) ) ), dCTF( fc == 16000 ) );
end

%% 5) cache per HRTF set
if numel( hH ) >= 2
    A = get_hrtf_equalization_filter( 'FF', hH(1) );
    B = get_hrtf_equalization_filter( 'FF', hH(2) );
    C = get_hrtf_equalization_filter( 'FF', hF );
    ok5 = ~isequal( A.left, B.left ) && ~isequal( A.left, C.left ) ...
          && strcmp( A.hrtf, hH(1).name ) && strcmp( B.hrtf, hH(2).name ) && strcmp( C.hrtf, hF.name );
    fprintf( '\n5) equalization filters cached per HRTF set (%s, %s, %s) -> %s\n', ...
             A.hrtf, B.hrtf, C.hrtf, il_pass( ok5 ) );
end

%% local functions
function n0 = il_onset( h )
% onset (samples) of each column: first sample within -20 dB of the maximum, after
% 10x upsampling (as in Brinkmann et al. 2019, without the low-pass filter)
hu = resample( h, 10, 1 );
n0 = zeros( 1, size( h, 2 ) );
for c = 1:size( h, 2 )
    a = abs( hu(:,c) );
    n0(c) = ( find( a >= max( a ) * 10^(-20/20), 1 ) - 1 ) / 10;
end
end

function t = il_pass( ok )
if ok, t = 'PASSED'; else, t = 'FAILED'; end
end
