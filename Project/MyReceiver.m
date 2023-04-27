function [OutPutBits,Feedback] =  MyReceiver(InputSamples, Parameters)

%% frequency synchronization
[FreqOffset]  =  FreqOffsetEstimate(exp(1j*2*pi*Parameters.freqOffset*(0:1:length(InputSamples)-1)) .* InputSamples, Parameters.MaxOffset, Parameters.NumOffsets, Parameters.h);
rxSamplesFreqOffsetCorr =  InputSamples .* exp(-1j*2*pi*FreqOffset*(0:1:length(InputSamples)-1));

%% time synchronization


%% Matched-filtered
rxSamplesMF = conv(rxSamplesFreqOffsetCorr,Parameters.h);
rxSamplesMFDelayCorr = rxSamplesMF(Parameters.delay:Parameters.sps:(Parameters.E/(log2(Parameters.M))-1)*Parameters.sps + Parameters.delay);
%% channel Estimator

%% demodulation
%% modulation
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