function [OutPutBits,Feedback] =  MyReceiver(InputSamples, Parameters)



%% frequency synchronization
[FreqOffset]  =  FreqOffsetEstimate(InputSamples, Parameters.MaxOffset, Parameters.NumOffsets, Parameters.h);
rxSigFreqCorr =  InputSamples .* exp(-1j*2*pi*FreqOffset*(0:1:length(InputSamples)-1));

%% time synchronization
%[tau] = SymbolTimingEstimate(rxSigFreqCorr, Parameters);



%% channel Estimator
if(Parameters.ChannelType == "AWGN")
    Ns = (Parameters.E+Parameters.numAppendBits)/(log2(Parameters.M));
    % Matched-filtered
    rxSamplesMF = conv(rxSigFreqCorr,Parameters.h);
    rxSamplesEq = rxSamplesMF(Parameters.delay:Parameters.sps:(Ns-1)*Parameters.sps + Parameters.delay);
else
    Ns = (Parameters.E+Parameters.numAppendBits)/(log2(Parameters.M)) + Parameters.Np;
    Np = (Parameters.Np);
    if(Parameters.PerfectChannelEst == "YS")
        % Matched-filtering
        rxSamplesMF = conv(rxSigFreqCorr .* (conj(Parameters.Channel)./(abs(Parameters.Channel))),Parameters.h);
        rxSamplesMFOpt = rxSamplesMF(Parameters.delay:Parameters.sps:(Ns-1)*Parameters.sps + Parameters.delay);
        %rxSamplesMFOpt = rxSamplesMF;
        rxSamplesEq    = rxSamplesMFOpt(Parameters.symLoc);
    else
        pilots = Parameters.pilots;
        pilotLoc = Parameters.pilotLoc;
        rxSamplesMF = conv(rxSigFreqCorr,Parameters.h);
        rxSamplesMFOpt = rxSamplesMF(Parameters.delay:Parameters.sps:(Ns-1)*Parameters.sps + Parameters.delay);
        %rxSamplesMFOpt = rxSamplesMF;
        gamma_p = conj(pilots.') .* rxSamplesMFOpt(pilotLoc);
        if(0)


            Gamma_p = fft(gamma_p);
            %Gamma = [Gamma_p(1:ceil((Np-1)/2)+1) zeros(1,Ns-Np) Gamma_p(floor((Np-1)/2) + 2:Np)];
            Gamma = [Gamma_p(1:2) zeros(1,Ns-Np) Gamma_p(3:Np)];
            gamma = ifft(Gamma) * Ns/Np;

            rxSamplesEq = rxSamplesMFOpt./gamma;
        elseif(1)
            W = zeros(Np, Ns);
            Phi_rr = zeros(Np,Np);
            for k = 1:Np
                for i = 1:Np
                    if k == i
                        Phi_rr(k,k) = besselj(0, 2*pi*Parameters.fd*abs(i-k)*(Parameters.NpS+1)*4/Parameters.fs) + 1/Parameters.SNR;
                    else
                        Phi_rr(k,i) = besselj(0, 2*pi*Parameters.fd*abs(i-k)*(Parameters.NpS+1)*4/Parameters.fs);
                    end
                end
            end
            for i = 1:Ns
                phi_rgamma = zeros(Np,1);
                for k = 1:Np
                    phi_rgamma(k) = besselj(0, 2*pi*Parameters.fd*abs(i-1+(k-1)*(Parameters.NpS+1))*4/Parameters.fs);
                end
                W(:,i) = Phi_rr^(-1)*phi_rgamma;
            end
            gamma = W'*gamma_p.';
            rxSamplesEq = rxSamplesMFOpt./gamma.';
        else
            
        end


        %         rxSamplesEqT = InputSamples .* (conj(Parameters.Channel)./(abs(Parameters.Channel).^2));
        %         rxSamplesMFT = rxSamplesEqT(Parameters.symLoc);

        rxSamplesEq = rxSamplesEq(Parameters.symLoc);
    end
    
end

%% demodulation
if Parameters.Modtype == "PSK"
    rxBitsllr = pskdemod(rxSamplesEq.',Parameters.M, pi/4 ,"gray","OutputType","approxllr","NoiseVariance", (1./Parameters.SNR).^2);
else
    rxBitsllr = qamdemod(rxSamplesEq.',Parameters.M,'UnitAveragePower', 1,"OutputType","approxllr","NoiseVariance", (1./Parameters.SNR).^2);
end
%% decoder
[OutPutBits, ~] = ChannelDecoder(rxBitsllr, Parameters);
%OutPutBits = decoderBits;
end