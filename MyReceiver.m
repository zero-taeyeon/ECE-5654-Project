function [OutPutBits,Feedback] =  MyReceiver(InputSamples, Parameters)

%% frequency synchronization
InputSamples = exp(1j*2*pi*Parameters.freqOffset*(0:1:length(InputSamples)-1)) .* InputSamples;
FreqOffset = FreqOffsetEstimate(InputSamples, Parameters.MaxOffset, Parameters.NumOffsets, Parameters.h);
FreqOffset = 0;
rxSamplesFreqOffsetCorr = InputSamples .* exp(-1j*2*pi*FreqOffset*(0:1:length(InputSamples)-1));

%% time synchronization


%% Matched-filtered
rxSamplesMF = conv(rxSamplesFreqOffsetCorr,Parameters.h);
rxSamplesMFDelayCorr = rxSamplesMF(Parameters.delay:Parameters.sps:((Parameters.E+round(Parameters.E/Parameters.N))/(log2(Parameters.M))-1)*Parameters.sps + Parameters.delay);
%% channel Estimator
rp = rxSamplesMFDelayCorr(1:21:length(rxSamplesMFDelayCorr));
gammap = rp/Parameters.pilot;
W = zeros(length(rp), 544);
Phi_rr = zeros(length(rp),length(rp));
for k = 1:length(rp)
    for i = 1:length(rp)
        if k == i
            Phi_rr(k,k) = besselj(0, 2*pi*Parameters.fd*abs(i-k)*(Parameters.N+1)) + 1/Parameters.SNR;
        else
            Phi_rr(k,i) = besselj(0, 2*pi*Parameters.fd*abs(i-k)*(Parameters.N+1));
        end
    end
end
for i = 1:544
    phi_rgamma = zeros(length(rp),1);
    for k = 1:length(rp)
        phi_rgamma(k) = besselj(0, 2*pi*Parameters.fd*abs(i-1+(k-1)*(Parameters.N+1)));
    end
    W(:,i) = Phi_rr^(-1)*phi_rgamma;
end

gamma = W'*gammap';
rxSamplesMFDelayCorr = rxSamplesMFDelayCorr./gamma';
rxSamplesMFDelayCorr(1:21:length(rxSamplesMFDelayCorr)) = [];
%% demodulation

if Parameters.Modtype == "PSK"
    %[b_cap] = MyDetectPSK(rxSamplesFreqOffsetCorr,Parameters.M);
    rxBitsllr = pskdemod(rxSamplesMFDelayCorr,Parameters.M,0,"gray","OutputType","llr","NoiseVariance", (1./Parameters.SNR).^2);
else
    rxBitsllr = qamdemod(rxSamplesMFDelayCorr,Parameters.M,"gray","OutputType","llr","NoiseVariance", (1./Parameters.SNR).^2);
    %[b_cap] = MyDetectQAM(rxSamplesFreqOffsetCorr,Parameters.M);
end
%% decoder
[decoderBits, Parameters] = ChannelDecoder(rxBitsllr, Parameters);
OutPutBits = decoderBits;
end