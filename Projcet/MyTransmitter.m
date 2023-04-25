function OutputSamples = MyTransmitter(b, Parameters)
%% modulation
if Parameters.Modtype == "PSK"
    [~,symbols] = MyPSK(b,Parameters.M); 
else
    [~,symbols] = MyQAM(b,Parameters.M); 
end
%% pulse shaping
sps = 4;            % samples for symbol
numTaps = 25;
rolloff   = 0.25;
% Assuming sample period is 0.1 micro sec, so sample rate is 10 Msps, 
% so *symbol* period is sps

if(Parameters.pulseShape == "SQAR")
    delay = sps;
else
    delay = sps*numTaps+1; %length(h);
end
[OutputSamples, h, Eg] = PulseShape(symbols, sps, numTaps, rolloff, Parameters.pulseShape);

end