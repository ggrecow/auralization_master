function varargout = ita_aurelio_control(varargin)
%ITA_AURELIO_CONTROL - Send Settings to Aurelio 2014 Frontend
%  This function sends Midi Sysex Commands to the Aurelio High Precision
%  Frontend (Swen Mueller, Immetro, Rio de Janiero, Brazil).
% 
%
%  Syntax: ita_aurelio_control(options)
%       'init' (false): set the frontend to last know values
%       'range' (0): inputRange as double, will round to nearest possible
%       'samplingRate' (ita_preferences), 32000, 44100 or 48000 multiplied by 1x, 2x, or 4x.
%       'coupling'
%
%   Reference page in Help browser
%        <a href="matlab:doc ita_robocontrol">doc ita_robocontrol</a>

% <ITA-Toolbox>
% This file is part of the application FrontendControl for the ITA-Toolbox. All rights reserved.
% You can find the license for this m-file in the application folder.
% </ITA-Toolbox>

% Author: Pascal Dietrich - pdi@akustik.rwth-aachen.de
% reworked: Jan Richter - jri@akustik.rwth-aachen.de
% reworked: Jonas Förster - jonas.foerster@rwth-aachen.de: commands for
% CMD42 added (11/2018): for a detailed description of these commands see
% the documentation for the CMF42 (CMF42 - product description.pdf)

%% persistent
% persistent oldSettings last_input last_AmpHighPower last_Amplifier last_AmpBridgeMode last_AmpLowImpedanceMode last_Amp26dBu last_AmpAC last_AmpMono
persistent settings presets presetNames currentPresetIndex;
% force_init = false;
if isempty(settings)  
    settings = getInitSettings();
end

if nargin == 0
    ita_aureliocontrol_gui();
end

%% init
sArgs = getArgsFromSettings(settings);

%% parse
[sArgs]   = ita_parse_arguments(sArgs,varargin);
ch_number = ita_angle2str( sArgs.channel - 1, 2); %zero indexing for channel numbering

if sArgs.getSettings
    varargout{1} = settings;
    varargout{2} = presetNames;
    varargout{3} = currentPresetIndex;
    if ~isempty(presets)
        presetChanged = ~isequal(presets{currentPresetIndex},settings);
    else
        presetChanged = 1;
    end
    varargout{4} = presetChanged;
    return;
end


if sArgs.savePreset
    presetIndex = length(presetNames) + 1;
    presets{presetIndex} = settings;
    presetNames{presetIndex} = sArgs.presetName;
    currentPresetIndex = presetIndex;
end


if sArgs.getPresets
    
end

if sArgs.reset
    clear all
    ita_aurelio_control('init');
    if nargout == 1
        varargout{1} = [];
    end
    return;
end

%% inverse parameters
if sArgs.Amp06dBu
    sArgs.Amp26dBu = false;
end

if sArgs.Amp26dBu == false
   sArgs.Amp06dBu = true; 
end

if sArgs.NoAmplifier
    sArgs.Amplifier = false;
end

if sArgs.Amplifier
    sArgs.NoAmplifier = false;
end

if sArgs.AmpLowPower
    sArgs.AmpHighPower = false;
end

if sArgs.AmpHighPower
    sArgs.AmpLowPower = false;
end

if sArgs.NoGroundLift
    sArgs.groundLift = false;
end

if sArgs.groundLift
    sArgs.NoGroundLift = false;
end

settings = setSettingsFromArgs(sArgs,settings);

if sArgs.setPreset
    % find the preset with the given name
    for index = 1:length(presetNames)
       if strcmp(presetNames{index},sArgs.presetName)
          currentPresetIndex = index;
          settings = presets{index};
       end
    end
end

%% write back to persistents
% last_input          = sArgs.input;
% last_AmpHighPower   = sArgs.AmpHighPower;
% last_Amplifier      = sArgs.Amplifier;
% last_AmpBridgeMode  = sArgs.AmpBridgeMode;
% last_AmpLowImpedanceMode = sArgs.AmpLowImpedanceMode;
% last_Amp26dBu       = sArgs.Amp26dBu;
% last_AmpAC          = sArgs.AmpAC;
% last_AmpMono        = sArgs.AmpMono;

%% INIT device
if sArgs.init
    % go thru all stages
    clear
    ita_aurelio_control('input','XLR','inputrange',6,'feed',0 , 'samplingRate',ita_preferences('samplingRate'));
    return
end

%% 00h: coupling and feed control
% if ~isempty( sArgs.feed )
for index = 1:length(sArgs.channel)
    rS = getCommandFromFeedSelect(sArgs,sArgs.channel(index),settings);

    par_value  = dec2hex( bin2dec ( ['0' rS.Wait rS.Lem28 rS.Phan rS.Feed rS.ICP rS.Glift rS.AC] ));
    send_sysex(rS.par_number, par_value, ita_angle2str(sArgs.channel(index)-1,2)); % send final sysex
end    
    
    
% end

%% 01h: input select -- routing control -- HUHU old values needed

par_number = '01';
switch lower(settings.mode)
    case 'norm' %normal mode
        Mode = '0000';
    case 'imp' %impedance measurement
        Mode = '0001';
        sArgs.feed = false;
        %             ChSwp = '1'; % pdi: oct 2012: for some reason
        %             this is not required anymore ???!!!
    case {'impref','iref'}
        Mode = '0010';
        sArgs.feed = false;
    case {'bncref','bref'}
        Mode = '0100';
        sArgs.input = 'gnd';
        sArgs.feed = false;
    case {'ampref','aref'}
        Mode = '0101';
        sArgs.input = 'gnd';
        sArgs.feed = false;
    case {'xlrref','xref','lineref'}
        Mode = '0011';
        sArgs.feed = false;
    case 'specialref'
        Mode = ['010' num2str(~ch_number) ];
        sArgs.input = 'gnd';
    otherwise
        error('argument not correct for ''input''')
end

if ~exist('ChSwp','var')
    ChSwp = '0'; %funny channel swapping for crazy people. take care, dude!
else
    ita_verbose_info('Careful, channel swapping is activated.',0)
end

% get input selection
IS = getCommandFromInputSelect(settings.ch(1).inputselect);
par_value = dec2hex( bin2dec( [IS ChSwp Mode ] ) );
send_sysex(par_number, par_value, ita_angle2str(0,2));

IS = getCommandFromInputSelect(settings.ch(2).inputselect);
par_value2 = dec2hex( bin2dec( [IS ChSwp Mode ] ) );
send_sysex(par_number, par_value2, ita_angle2str(1,2)); %send to device

%% 02h: range

for index = 1:length(sArgs.channel)
    sArgs.inputrange = min( max(settings.ch(sArgs.channel(index)).inputrange,-34) , 56);
    par_number = '02';
    par_value  = round((- sArgs.inputrange + 56)/10); %round to nearest possible
%         for idx = 1:numel(sArgs.channel)
%             disp(idx)
%             56 - (par_value * 10)
%             settings.ch(sArgs.channel(idx)).inputrange = 56 - (par_value * 10);
%         encd
    if sArgs.securityMode
        par_value  = par_value + 40;
    end
    par_value = ita_angle2str(par_value,2);
    send_sysex(par_number, par_value, ita_angle2str(sArgs.channel(index)-1,2)); %send to device
end

%% 03h: sampling rate

par_number = '03';

if isnatural(sArgs.samplingRate / 48000)
    modifier  = (sArgs.samplingRate / 48000);
    base_rate = 2;
elseif isnatural(sArgs.samplingRate / 44100)
    modifier  = sArgs.samplingRate / 44100;
    base_rate = 1;
elseif isnatural(sArgs.samplingRate / 32000)
    modifier  = sArgs.samplingRate / 32000;
    base_rate = 0;
else
    error('sampling rate not supported.')
end
switch modifier
    case 1
        modifier = 0;
    case 2
        modifier = 1;
    case 4
        modifier = 2;
    otherwise
        error('sampling rate not supported.')
end
par_value  =  dec2hex(modifier*4 + base_rate);
send_sysex(par_number, par_value, []); % send final sysex

ita_preferences('samplingRate',sArgs.samplingRate); % also initializes playrec

%% 04h: Input Section Digital I/O Control
par_number = '04';
Dig = num2str(settings.Dig);
Pro = num2str(settings.Pro);
Sync1 = num2str(settings.Sync1); % TODO: what do these two values mean?
Sync0 = num2str(settings.Sync0);
DIn = settings.DIn;

par_value = dec2hex( bin2dec( [Sync1 Sync0 Pro Dig DIn] ) );

send_sysex(par_number, par_value, []); %send to device


%% 05h: analog output control - Parameter 05
HiPow   = num2str(settings.AmpHighPower);
AmpOn   = num2str(settings.Amplifier);
Bridge  = num2str(settings.ampBridgeMode);
LoImp   = num2str(settings.ampLowImpedanceMode);
dBu26   = num2str(settings.Amp26dBu);
AC      = num2str(settings.ampAC);
mono    = num2str(settings.ampMono);

par_number = '05';
par_value  = dec2hex(bin2dec(['0' HiPow AmpOn Bridge LoImp dBu26 AC mono])); % 0dB attenuation
send_sysex(par_number, par_value, ch_number); %send to device

if settings.Amp26dBu
    settings.amp_gain = 0;
else
    settings.amp_gain = -20;
end

%% 06h: output volume - Parameter 06
% if ~isempty( sArgs.outputvolume )
%     par_number = '06';
%     par_value  = '7F'; % 0dB attenuation
%     send_sysex(par_number, '0', []); %send to device
% end

%% 07h: Output Section Digital I/O Control
par_number = '07';
Ana = num2str(settings.Ana);
OutOn = num2str(settings.outOn);
Slow = num2str(settings.Slow);
FPCon = num2str(settings.FPCon);
Auto = num2str(settings.Auto);
Din1 = num2str(settings.Din1); %TODO: change to single array
Din0 = num2str(settings.Din0);

par_value = dec2hex( bin2dec( [Ana OutOn Slow FPCon Auto Din1 Din0] ) );
send_sysex(par_number, par_value, []); %send to device

%% 0Ah: set "AD Level Autosend"-properties of the Aurelio
% Function which tells Aurelio to periodically send the 
% AD/DA levels without request.

par_number = '0A';

if settings.Autosend >= 0 && settings.Autosend <= 1270
autosend = settings.Autosend;
autosendBin = dec2bin(autosend /10);
% h00 <- disable autosend property; otherwise: send AD/DA levels
% periodically after a defined time (8 bit binary number x 10 ms)
par_value  = dec2hex( bin2dec ( autosendBin ));

send_sysex(par_number, par_value, []); %send to device
end

%% 0Eh: programmable gain ampflifier (PGA) sensitivity

par_number = '0E';

PGAsens = settings.PGAsens; %in dB
if PGAsens >=0 && PGAsens <= 31.5
    PGAsensDec = PGAsens * 2;
    PGAsensBin = dec2bin(PGAsensDec); % round
    par_value = dec2hex( bin2dec ( ['00' PGAsensBin] ));
    send_sysex(par_number, par_value, ch_number); %send to device
else
    ita_verbose_info('Value for the PGA sensitivity is not valid!',0)
end

%% 11h: Power amplifier input control (v7 and above)
par_number = '11';

Use17 = num2str(settings.Use17); % compatibility mode to older versions
FourCH = num2str(settings.fourCH); % on: RME channles 1-4 used, off: RME channels 1-2 used
min10dB = num2str(settings.min10dB); % activate the -10dB pad in the input stage of the power amp
AC = num2str(settings.AC); % activate first order highpass in the input stage 
inputSelect = settings.inputSelect; % 

switch inputSelect
    case 'DAC Ch1/2'
        Mix = '0';
        Ana = '0';
    case 'ANA D15'
        Mix = '0';
        Ana = '1';
    case 'DAC Ch1/2 + 3/4'
        Mix = '1'; % default 2-ch mode
        Ana = '0';
    case 'ANA + DAC Ch 3/4'
        Mix = '1';
        Ana = '1';
end

par_value = dec2hex( bin2dec ( ['0' Use17 '0' FourCH min10dB AC Mix Ana] ));

send_sysex(par_number, par_value, []); %send to device

end

%% **********************  send sysex  ************************************
function send_sysex(par_number, par_value, ch_number)

% disp([par_number par_value ch_number])
% 
if ~isempty(ch_number) && numel(str2num(ch_number)) >= 2
    ch_number = '7F';
end
%build complete sysex
sys_hex = {};
sys_hex{numel(sys_hex)+1} = par_number; % hex
sys_hex{numel(sys_hex)+1} = par_value;  % hex
if ~isempty(ch_number)
    sys_hex{numel(sys_hex)+1} = ch_number;  % hex
end

for idx = 1:numel(sys_hex)
    sys_dec(idx) = hex2dec( sys_hex{idx} );
end

% generate checksum
checksum = sum(sys_dec);
checksum = bin2dec(num2str(mod(str2double(dec2bin(checksum)), 10000000)));
complete_sysex = [sys_dec checksum];

%send sysex
sysex = [hex2dec('F0') hex2dec('70') complete_sysex hex2dec('F7')]; %pre- post- ampel
ita_midi(sysex,ita_preferences('out_midi_DeviceID'));

end


%% get settings

function settings = getInitSettings()

    settings.AmpHighPower = false;
    settings.Amplifier = true;
    settings.ampBridgeMode = false;
    settings.ampLowImpedanceMode = false;
    settings.Amp26dBu = false;
    settings.ampAC = true;
    settings.ampMono = false;
    settings.samplingRate = ita_preferences('samplingRate');
    settings.mode = 'norm';
    settings.groundLift = false;
    settings.inputCouplingAC = true;
    settings.ch(1).inputrange = 6;
    settings.ch(2).inputrange = 6;
    settings.ch(1).inputfeed = 'none';
    settings.ch(2).inputfeed = 'none';
    settings.ch(1).inputselect = 'xlr';
    settings.ch(2).inputselect = 'xlr';
    %------------------------ new commands (11/18)-------------------------
    settings.Ana = false;
    settings.outOn = false;
    settings.Slow = false;
    settings.FPCon = false;
    settings.Auto = true;
    settings.Din1 = true;
    settings.Din0 = false;
    settings.Dig = false;
    settings.Pro = false;
    settings.Sync1 = false;
    settings.Sync0 = false;
    settings.DIn = '00';
    settings.PGAsens = 0;
    % ---------------- CMF42 ----------------------------------------------
    settings.Use17 = false; % false = compatibility mode to CMF22
    settings.fourCH = false;
    settings.min10dB = true;
    settings.AC = true;
    settings.inputSelect = 'DAC Ch1/2 + 3/4';
    settings.Autosend = 0;
%     settings.amp_gain = -20;

end

function sArgs = getArgsFromSettings(settings)

sArgs     = struct('input',[],'feed',[],'inputrange',[],'channel',1:2,'groundLift',settings.groundLift,'inputCouplingAC',settings.inputCouplingAC,...
    'mode',settings.mode,'securityMode',false,'init',false,'samplingRate',settings.samplingRate,'outputvolume',[],...
    'AmpHighPower',settings.AmpHighPower,'AmpLowPower',false,'Amplifier',settings.Amplifier,'NoAmplifier',false,'AmpBridgeMode',settings.ampBridgeMode,'AmpLowImpedanceMode',settings.ampLowImpedanceMode,...
    'Amp26dBu',settings.Amp26dBu,'Amp06dBu',false , 'NoGroundLift', false, 'AmpAC', settings.ampAC, 'AmpMono',...
    settings.ampMono, 'Ana', settings.Ana,'outOn', settings.outOn, 'Slow', settings.Slow,'FPCon', settings.FPCon,...
    'Auto', settings.Auto, 'Din1', settings.Din1, 'Din0', settings.Din0, 'Dig', settings.Dig, 'Pro', settings.Pro,...
    'Sync1', settings.Sync1, 'Sync0', settings.Sync0, 'DIn', settings.DIn, 'PGAsens', settings.PGAsens,...
    'Use17', settings.Use17, 'fourCH', settings.fourCH, 'min10dB', settings.min10dB, 'AC', settings.AC,...
    'inputSelect', settings.inputSelect,'Autosend', settings.Autosend,...
    'getSettings',false,'reset',false,'setPreset',false,'getPresets',false,'savePreset',false,'presetName','');

end

function settings = setSettingsFromArgs(sArgs,settings)

    settings.AmpHighPower = sArgs.AmpHighPower;
    settings.Amplifier = sArgs.Amplifier;
    settings.ampBridgeMode = sArgs.AmpBridgeMode;
    settings.ampLowImpedanceMode = sArgs.AmpLowImpedanceMode;
    settings.Amp26dBu = sArgs.Amp26dBu;
    settings.ampAC = sArgs.AmpAC;
    settings.ampMono = sArgs.AmpMono;
    settings.samplingRate = sArgs.samplingRate;
    settings.mode = sArgs.mode;
    settings.groundLift = sArgs.groundLift;
    settings.inputCouplingAC = sArgs.inputCouplingAC;
    settings.Ana = sArgs.Ana;
    settings.outOn = sArgs.outOn;
    settings.Slow = sArgs.Slow;
    settings.FPCon = sArgs.FPCon;
    settings.Auto = sArgs.Auto;
    settings.Din1 = sArgs.Din1;
    settings.Din0 = sArgs.Din0;
    settings.Dig = sArgs.Dig;
    settings.Pro = sArgs.Pro;
    settings.Sync1 = sArgs.Sync1;
    settings.Sync0 = sArgs.Sync0;
    settings.DIn = sArgs.DIn;
    settings.PGAsens = sArgs.PGAsens;
    settings.Use17 = sArgs.Use17;
    settings.fourCH = sArgs.fourCH;
    settings.min10dB = sArgs.min10dB;
    settings.AC = sArgs.AC;
    settings.inputSelect = sArgs.inputSelect;
    settings.Autosend = sArgs.Autosend;
    for index = 1:length(sArgs.channel)
        if ~isempty(sArgs.inputrange)
            settings.ch(sArgs.channel(index)).inputrange = sArgs.inputrange;
        end
        if ~isempty(sArgs.feed) 
            settings.ch(sArgs.channel(index)).inputfeed = sArgs.feed;
        end
        if ~isempty(sArgs.input)
            settings.ch(sArgs.channel(index)).inputselect = sArgs.input;
        end
    end
end

% input selection
function IS = getCommandFromInputSelect(in)
    
    switch lower(in)
        case 'xlr'
            IS = '11';
        case 'lemo'
            IS = '01';
            sArgs.feed = 'pol';
        case 'gnd'
            IS = '00';
        case 'bnc'
            IS = '10';
        otherwise
            %JRI error
            disp(['input select unknown'])
    end

end

function returnStruct = getCommandFromFeedSelect(sArgs,channel,settings)
    par_number = '00';
    Wait  = '0'; %wait for relays to switch later
    Lem28 = '0'; %switch 14 to 28Volts, Pin7 is then grounded
    Phan  = '0';
    Feed  = '0';
    ICP   = '0';
    Glift = num2str(settings.groundLift);
    AC    = num2str(settings.inputCouplingAC);
    
    switch lower(settings.ch(channel).inputfeed)
        case 'pha'
            Phan = '1';
            AC   = '1'; %block DC from preamp inputs
        case 'pol'
            Feed = '1';
            AC   = '1'; %block DC from preamp inputs
        case {'icp','iepe'}
            ICP  = '1';
            Feed = '1';
            AC   = '1';
        case 'p+p'
            Phan = '1';
            Feed = '1';
            AC   = '1';
        case 'all'
            Phan = '1';
            Feed = '1';
            ICP  = '1';
            AC   = '1';
        case 'ccx'
            Phan = '0';
            Feed = '0';
            ICP  = '0';
            AC   = '0';
            Glift = '1';
        case {0 ,'none','off'}
            %
        otherwise
            error('feed wrong')
    end
    
    returnStruct.par_number = par_number;
    returnStruct.Wait = Wait;
    returnStruct.Lem28 = Lem28;
    returnStruct.Phan = Phan;
    returnStruct.Feed = Feed;
    returnStruct.ICP = ICP;
    returnStruct.Glift = Glift;
    returnStruct.AC = AC;
end

