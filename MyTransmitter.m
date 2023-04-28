function [OutputSamples, Parameters] = MyTransmitter(b, Parameters)
%% forward error correction code

[encoderBits, Parameters] = ChannelEncoder(b, Parameters);

%% adding pilot symbol

%% modulation
if Parameters.Modtype == "PSK"
    %[~,symbols] = MyPSK(encoderBits,Parameters.M); 
    unpadded_symbol = pskmod(encoderBits, Parameters.M);
else
    [~,unpadded_symbol] = MyQAM(encoderBits,Parameters.M); 
end
Parameters.pilot = pskmod(zeros(log2(Parameters.M),1), Parameters.M, InputType='bit');
N = 20;
unpadded_symbol = unpadded_symbol';
symbols = [];
for i = 1:N:518
    if 518-i < N
        symbols = [symbols Parameters.pilot unpadded_symbol(i:end)];
    else
        symbols = [symbols Parameters.pilot unpadded_symbol(i:i+N-1)];
    end
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