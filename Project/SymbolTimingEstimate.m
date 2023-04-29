function [tau]  =  SymbolTimingEstimate(InputSignal, Parameters)
sps = Parameters.sps;
N = Parameters.Ns;
delay = length(Parameters.h);
sampleT = delay:1:(N-1)*sps + delay + 3;
rxSamplesABS = (abs(InputSignal(sampleT))).^2 .*  exp(-1j*2*pi*(0:1:length(sampleT)-1)/(sps));
tau  = (-(sps/(2*pi))*unwrap(angle(sum(rxSamplesABS))));
if(tau < 0)
    tau = mod(tau,4);
end
end