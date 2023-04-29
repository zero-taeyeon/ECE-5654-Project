function [OutPutBits,Feedback] =  MyReceiver(InputSamples, Parameters)

Ns = (Parameters.E+Parameters.numAppendBits)/(log2(Parameters.M)) + Parameters.Np;
Np = (Parameters.Np);

%% frequency synchronization
% [FreqOffset]  =  FreqOffsetEstimate(exp(1j*2*pi*Parameters.freqOffset*(0:1:length(InputSamples)-1)) .* InputSamples, Parameters.MaxOffset, Parameters.NumOffsets, Parameters.h);
% rxSamplesFreqOffsetCorr =  InputSamples .* exp(-1j*2*pi*0*(0:1:length(InputSamples)-1));

%% time synchronization

%% channel Estimator
if(Parameters.ChannelType == "AWGN")
    % Matched-filtered
    rxSamplesMF = conv(rxSamplesFreqOffsetCorr,Parameters.h);
    rxSamplesMFOpt = rxSamplesMF(Parameters.delay:Parameters.sps:(Ns-1)*Parameters.sps + Parameters.delay);
else
    if(Parameters.PerfectChannelEst == "YS")
        % Matched-filtering
        rxSamplesMF = conv(rxSamplesFreqOffsetCorr .* (conj(Parameters.Channel)./(abs(Parameters.Channel))),Parameters.h);
        rxSamplesMFOpt = rxSamplesMF(Parameters.delay:Parameters.sps:(Ns-1)*Parameters.sps + Parameters.delay);
    else
        if(0)
            % Matched-filtering
            rxSamplesMF = conv(InputSamples,Parameters.h);
            rxSamplesMFOpt = rxSamplesMF;
            %rxSamplesMFOpt = rxSamplesMF(Parameters.delay:Parameters.sps:(Ns-1)*Parameters.sps + Parameters.delay);
            % FFT Channel Estimation
            rxSamplesPilots = conj(Parameters.pilots.') .* rxSamplesMFOpt(Parameters.pilotLoc);
            pilotChfft = fft(rxSamplesPilots);
            pilotCh    = [pilotChfft(1:floor((Np-1)/2 + 1)) zeros(1, Ns-Np) pilotChfft(ceil((Np-1)/2 ) + 1 : Np)];
            chEst      = ifft(pilotCh) * Ns/Np ;
            chEst = interp1(Parameters.pilotLoc,rxSamplesPilots,setdiff(1:Ns,Parameters.pilotLoc),"linear","extrap");
            H = Parameters.Channel;
            rxSamplesMFOpt =  rxSamplesMFOpt .* conj(chEst)./(abs(chEst));
            rxSamplesMFOpt = rxSamplesMFOpt(Parameters.symLoc);
        else

            % channel estimation
            p = Parameters.pilots;
            P = zeros(Ns,1);
            P(Parameters.pilotLoc) = p;
            % Apply pulse shaping
            pilots_PulseShaped = PulseShape(P,Parameters.sps,Parameters.numTaps,Parameters.rolloff,Parameters.pulseShape);
            rxSamplesPilots = conj(Parameters.pilots.') .* InputSamples(Parameters.pilotLoc);
            h_LS = zeros(1,length(pilots_PulseShaped));
            h_LS(Parameters.pilotLoc) = InputSamples(Parameters.pilotLoc) .* (conj(InputSamples(Parameters.pilotLoc)));

            rxSamplesEq = ifft(fft(InputSamples).*conj(fft(h_LS))./abs(fft(h_LS)));
            rxSamplesMFOpt = rxSamplesEq(Parameters.symLoc);
           
            
            
            
            
            
            %h_LS = rxSamplesFreqOffsetCorr .* (conj(pilots_PulseShaped)./abs(pilots_PulseShaped).^2);
            % pilot_peak = (Parameters.delay-1)/2+1:Parameters.NpS*Parameters.sps:Parameters.NpS*Parameters.Np*Parameters.sps+(Parameters.delay-1)/2;
            %h_LS = zeros(1,length(pilots_PulseShaped));
%             gamma_p = rxSamplesFreqOffsetCorr(pilot_peak) ...
%                 .* conj(pilots_PulseShaped(pilot_peak))./abs(pilots_PulseShaped(pilot_peak)).^2;
%             














            %             h_LS = rxSamplesFreqOffsetCorr .* conj(pilots_PulseShaped)./abs(pilots_PulseShaped()).^2;
            %             rxSamplesEq = ifft(fft(rxSamplesFreqOffsetCorr).*conj(fft(h_LS))./abs(fft(h_LS)));
            %             pilot_peak = (Parameters.delay-1)/2+1:Parameters.NpS*Parameters.sps:Parameters.NpS*Parameters.Np*Parameters.sps;
            %             h_LS = zeros(1,length(pilots_PulseShaped));
            %             h_LS(pilot_peak) = rxSamplesFreqOffsetCorr(pilot_peak) ...
            %                 .* conj(pilots_PulseShaped(pilot_peak))./abs(pilots_PulseShaped(pilot_peak)).^2;
            %             %rxSamplesEq = ifft(fft(rxSamplesFreqOffsetCorr).*conj(fft(h_LS))./abs(fft(h_LS)));
%             H = interp1(pilot_peak,gamma_p(pilot_peak),setdiff(1:length(pilots_PulseShaped),pilot_peak),"linear","extrap");
            %             %rxSamplesMF = conv(rxSamplesEq,Parameters.h);
            %             rxSamplesMFOpt = rxSamplesEq((Parameters.delay-1)/2+1:Parameters.sps:(Ns)*Parameters.sps + (Parameters.delay-1)/2);
            %
%             H_cap = zeros(1,length(pilots_PulseShaped));
%             H_cap(pilot_peak) = gamma_p(pilot_peak);
%             H_cap(setdiff(1:length(pilots_PulseShaped),pilot_peak)) = H;


















            %H;
            %             rxSamplesMF = conv(rxSamplesEq,Parameters.h);
            %             rxSamplesMFOpt = rxSamplesMF(Parameters.delay:Parameters.sps:(Ns-1)*Parameters.sps + Parameters.delay);
            %            rxSamplesMFOpt = rxSamplesMFOpt(Parameters.symLoc);
            %             rp = rxSamplesMFOpt(Parameters.pilotVector).';
            %             p = Parameters.pilots;
            %             rxSamplesPilots = rp./p;
            %             phi_RR = zeros(Np,Np);
            %             for i = 1:1:Np
            %                 for k = 1:1:Np
            %                     %phi_RR(i,k) = besselj(0, 2*pi*Parameters.fd*abs(i-k)*Parameters.N./(Parameters.fs/4));
            %                     phi_RR(i,k) = besselj(0, 2*pi*Parameters.fd*abs(i-k)*(Parameters.N + 1) );
            %                 end
            %             end
            %             for i=1:1:Np
            %                 phi_RR(i,i) = phi_RR(i,i) + (1/Parameters.SNR);
            %             end
            %
            %             W = zeros(Np,Ns);
            %             k = (1:1:Np).';
            %             for i = 1:1:Ns
            %             %W(:,i) = inv(phi_RR) * besselj(0, 2*pi*Parameters.fd*abs(i-k*Parameters.N)./(Parameters.fs/4));
            %             W(:,i) = (phi_RR^(-1)) * besselj(0, 2*pi*Parameters.fd*abs(i-1+k*(Parameters.N+1)));
            %             end
            %
            %             chEst = (W' * rxSamplesPilots).';
            %             rxSamplesMFOpt =  rxSamplesMFOpt .* conj(chEst)./(abs(chEst));
        end
        %

    end
end

%% demodulation
if Parameters.Modtype == "PSK"
    rxBitsllr = pskdemod(rxSamplesMFOpt.',Parameters.M, 0 ,"gray","OutputType","llr","NoiseVariance", (1./Parameters.SNR).^2);
else
    rxBitsllr = qamdemod(rxSamplesMFOpt.',Parameters.M,'UnitAveragePower', 1,"OutputType","llr","NoiseVariance", (1./Parameters.SNR).^2);
end
%% decoder
[OutPutBits, ~] = ChannelDecoder(rxBitsllr, Parameters);
%OutPutBits = decoderBits;
end