function [OutputSamples, Parameters] = MyTransmitter(b, Parameters)

%% Forward error correction code
[encoderBits, Parameters] = ChannelEncoder(b, Parameters);

%% Modulation
if Parameters.Modtype == "PSK"
    [encoderBitsApp,symbols] = MyPSK(encoderBits,Parameters.M);
else
    [encoderBitsApp,symbols] = MyQAM(encoderBits,Parameters.M);
end

% Appended Bits
Parameters.numAppendBits = length(encoderBitsApp)-length(encoderBits);

if(Parameters.PerfectChannelEst == "YS")
    txSymbols = symbols;
else
    %% Pilot Addition
    % Pilot Parameters
    Parameters.NpS  = 11;
    Parameters.pilotLoc = 1:Parameters.NpS:length(symbols);
    Parameters.Np = length(Parameters.pilotLoc);
    [~,Parameters.pilots] = MyPSK(randi([0 1],Parameters.Np*2,1),4);

    txSymbols = zeros(length(symbols) + Parameters.Np,1);
    txSymbols(Parameters.pilotLoc) = Parameters.pilots;
    Parameters.symLoc = setdiff(1:length(symbols)+Parameters.Np,Parameters.pilotLoc);
    txSymbols(Parameters.symLoc) = symbols;
end
%% Pulse shaping
% Assuming sample period is 0.1 micro sec, so sample rate is 10 Msps,
% so *symbol* period is 0.4 micro secs consisting of four samples

if(Parameters.pulseShape == "SQAR")
    Parameters.delay = Parameters.sps;
else
    Parameters.delay = Parameters.sps*Parameters.numTaps+1; %length(h);
end

[OutputSamples, h, ~] = PulseShape(txSymbols, Parameters.sps, Parameters.numTaps, Parameters.rolloff, Parameters.pulseShape);

Parameters.h = h;       % pulse taps

end