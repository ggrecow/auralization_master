function [ OutputSignals ] = ita_hoa_decode( Bformat, LoudspeakerPos, varargin )
%ITA_DECODEAMBISONICS Decodes a BFormat Signal in ANC with SN3D
%normalization.
%   The term Bformat corresponds to the ambisonics
%   transimission format and is not bound to first order ambisonics in the
%   classic definition. Bformat can be either numeric or itaAudio.
%   LoudspeakerPos has to be itaCoordinates with nPoints.
%   Truncation order is determined by the number of channels of the input
%   format

%  BFormat<nmax,LS>

opts.decoding='remax'; % Decoding strategy (remax,inphase,plane)
opts.distanceloss=false; % Indicates to compensate the loudspeaker distance for a fully calibrated reproduction, not adapting to the listener position yet
opts.vrlab=false; % Mixes the virtual up and down loudspeaker to the upper and lower ring
opts.returnDecoder=false; %if you want to return the decoder the bformat should ones((N+1)^2,1);
opts = ita_parse_arguments(opts,varargin);

% Initializing further parameters
if isa(Bformat, 'itaAudio')
    nmax=max(Bformat.nChannels);
else
    nmax=size(Bformat,2);
end
N=floor(sqrt(nmax)-1);

% Weighting for Inphase/ReMax
g_re=ones(nmax,1); %init weighting factors for ReMax
g_in=g_re; %init weighting factors for InPhase


% ReMax Decoding
if(sum(strcmp(lower(opts.decoding),{'remax' 'both'})))
    %             Calculates the B-Format channel weights for ReMax-Decoding
    %             see J.Daniel Thesis p.312 A.67/A.68
    weightsReMax=zeros(nmax,1); % init weights
    % g_m=P_m(r_E) with n=0 for P
    % r_E is largest root of P_(order+1)
    syms x;
    
    f=1/2^(N+1)/factorial(N+1)*diff(((x^2-1)^(N+1)),(N+1));%Legendre Polynom n=0 m=order+1
    maxroot = double(max(vpasolve(f))); %find maximum root(Nullstelle)
    for k=1:nmax
        leggie=legendre(ita_sph_linear2degreeorder(k),abs(maxroot)); % g_m=P_m(r_E)
        weightsReMax(k)=leggie(1);% pick n=0
    end
    g_re=weightsReMax;
end

% Inphase Decoding
if(sum(strcmp(lower(opts.decoding),{'inphase' 'both'})))
    %             Calculates the B-format channel weights for InPhase-Decoding
    weightsInPhase=zeros(nmax,1);
    
    %Dissertation J.Daniel p.314, 'preserve Energy' in 0. order
    %g_0=sqrt[N*(2*M+1)/(M+1)^2]
    weightsInPhase(1)=sqrt(N*(2*N+1)/(N+1)^2);
    
    for k=2:nmax
        m=ita_sph_linear2degreeorder(k);
        % Dissertation J.Daniel p.314
        % g_m=M!*(M+1)!/[(M+m+1)!*(M-n)!]
        weightsInPhase(k)=factorial(N)*factorial(N+1)/factorial(m+N+1)/factorial(N-m);
    end
    g_in=weightsInPhase;
end

weights=g_in.*g_re; % merge weighting factors

%% Applying weighting to BFormat
if isa(Bformat,'itaAudio')
    for k=1:Bformat.nChannels
        Bformat.time(:,k)=Bformat.time(:,k).*weights(k);
    end
else
    for k=1:numel(weights)
        Bformat(:,k)=weights(k).*Bformat(:,k);
    end
end

%% SH and inversion of loudspeaker set-up
Y = ita_sph_base(LoudspeakerPos, N, 'real'); % generate basefunctions
Yinv=pinv(Y); % calculate Pseudoinverse, moore penrose, svd

if isa(Bformat,'itaAudio')
    OutputSignals=itaAudio(1,LoudspeakerPos.nPoints);
    for k=1:LoudspeakerPos.nPoints
        temp=itaAudio(nmax,1);
        for l=1:nmax
            temp(l)=Bformat.ch(l)*Yinv(l,k);
        end
        OutputSignals(k)=sum(temp);
        clear temp;
        
    end
    if opts.distanceloss || opts.vrlab
        for k=1:LoudspeakerPos.nPoints
            if opts.vrlab
                if k>4 && k<9
                    OutputSignals(k)=OutputSignals(k)+OutputSignals(14)*0.5;
                end
                if k>8 && k<13
                    OutputSignals(k)=OutputSignals(k)+OutputSignals(13)*0.5;
                end
            end
            if opts.distanceloss
                OutputSignals(k)=OutputSignals(k)*LoudspeakerPos.r(k);
            end
        end
    end
    if opts.vrlab
        OutputSignals(13:14)=[];
    end
    OutputSignals=ita_merge(OutputSignals(:));
    OutputSignals.channelCoordinates=LoudspeakerPos;
    
    
else
    OutputSignals=(Bformat*Yinv);
    if opts.vrlab %Distribute virtual loudspeaker up to 4 loudspeaker, (L1)^2+(L2)^2+(L3)^2+(L4)^2=1
        for k=5:8
            OutputSignals(:,k)=OutputSignals(:,k)+OutputSignals(:,14)*0.5;
        end
        for k=9:12
            OutputSignals(:,k)=OutputSignals(:,k)+OutputSignals(:,13)*0.5;
        end
    end
    if opts.distanceloss
        for k=1:numel(OutputSignals)
            OutputSignals(k)=OutputSignals(k)*LoudspeakerPos.r(k);
        end
    end
    if opts.vrlab
        OutputSignals(:,13:14)=[];
    end
    if ~opts.returnDecoder
        sum(OutputSignals,1);
    end
end


end
