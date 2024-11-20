function [ data, samplerate, metadata ] = dfRavenBinauralVA( alpha, beta, rpf )

    viewVec = itaCoordinates( 1 );
    viewVec.cart = [ sind(beta)*cosd(alpha) sind(beta)*sind(alpha) (-1)*cosd(beta)];

    upVec = viewVec;
    upVec.elevation= 90-viewVec.elevation;
    upVec.azimuth = 180+viewVec.azimuth;

    ravenRecViewVector = viewVec.cart([1 3 2]).*[1 1 -1];
    ravenRecUpVector = upVec.cart([1 3 2]).*[1 1 -1];

    rpf.setReceiverViewVectors(ravenRecViewVector);
    rpf.setReceiverUpVectors(ravenRecUpVector);
    % disp([ 'Alpha: ' sprintf('%0.2f ',alpha) ' | Beta: ' sprintf('%0.2f ',beta) ' |  View: ' sprintf('%0.2f ',ravenRecViewVector) ' | Up: ' sprintf('%0.2f ',ravenRecUpVector) ])
    % data = zeros(2,256);
    % samplerate = 44100;
    % metadata = [];

    rpf.run();
    brir = rpf.getBinauralImpulseResponseImageSourcesItaAudio;  % only early reflections are requested, use rpf.getBinauralImpulseResponseItaAudio if you want to get the full impulse response (and adjust the filter length in your rpf accordingly)

    n_residual = mod( brir.nSamples, 4 );
    data = [ brir.timeData', zeros( brir.nChannels, n_residual ) ];

    samplerate = brir.samplingRate;
    metadata = [];
end
