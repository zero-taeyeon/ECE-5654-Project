function [OutputSamples, Parameters] = MyTransmitter(b, Parameters)
%% forward error correction code

[encoderBits, Parameters] = ChannelEncoder(b, Parameters);

%% modulation
if Parameters.Modtype == "PSK"
    [~,symbols] = MyPSK(encoderBits,Parameters.M); 
else
    [~,symbols] = MyQAM(encoderBits,Parameters.M); 
end
%% pulse shaping

% Assuming sample period is 0.1 micro sec, so sample rate is 10 Msps, 
% so *symbol* period is sps

if(Parameters.pulseShape == "SQAR")
    Parameters.delay = Parameters.sps;
else
    Parameters.delay = Parameters.sps*Parameters.numTaps+1; %length(h);
end
[OutputSamples, h, Eg] = PulseShape(symbols, Parameters.sps, Parameters.numTaps, Parameters.rolloff, Parameters.pulseShape);

Parameters.h = h;

end