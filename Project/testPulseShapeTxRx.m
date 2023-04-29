clc;
clear;
close all;

%% pulse shape parameters
sps = 4;            % samples for symbol
numTaps = 25;
rolloff   = .25;
pulseShape = 'SqRa';
% Assuming sample rate is 1 Hz, so sample period is 1 sec,
% so *symbol* period is sps, Symbol Time = 8 sec => BW = 1/8
fs = 1;
f0 = 1/(sps*fs);

if(pulseShape == "SQAR")
    delay = sps;
else
    delay = sps*numTaps+1;%length(h);
end


%%% Parameters that can be changed by the student  %%%%
Parameters.ChannelType = 'AWGN'; % 'AWGN', 'RAYL', or 'RICE'
Parameters.fd = 10;              % max Doppler (10 or 100) only used for fading
Parameters.K = 10;               % Ricean K-factor; only used in Ricean fading
Parameters.PayloadSize = 100;   % this is completely up to the student
% This is the number of info bits per packet.
% The total packet size depends on the
% overhead included and coding rate.

Parameters.PerfectChannelEst = 'YS';% 'YS' - puts channel info into
% Parameters.Channel
% 'NO' - channel estimation must be done
% separately

%%% Parameters that SHOULD NOT bbe changed by the student  %%%%
Parameters.NumSinusoids = 100;      % used by Channel function
Parameters.fs = 1e7;                % used by Channel function


%% modulation Params
% generate bits
M = 2;
EbNodB = 0:20 ;
R = 1; % uncoded BPSK (1 bit/symbol)
SNR_dB = EbNodB + 10*log10(log2(M));
SNR = 10.^(SNR_dB/10);
EbNo = 10.^(EbNodB/10);
sigma = sqrt((10.^(-SNR_dB/10)));
N = 10; % number of bits of message per block
numBlocks = 200;
Nerrs = zeros(size(EbNodB));
offset = 0;

for j = 1:1:numBlocks

    txBits = randi([0 1],N,1);
    [txBits,symbols] = MyPSK(txBits,M);
    % pulse shaping
    [symbols_pulseShaped, h, Eg] = PulseShape(symbols, sps, numTaps, f0, rolloff, pulseShape);
    %symbols_pulseShaped = symbols_pulseShaped()
    %symbols_pulseShaped((sps*numTaps)/2+1+offset:sps:(length(symbols))*sps + (sps*numTaps)/2 +offset)
    for i = 1:1:length(EbNo)
        Parameters.SNR = SNR(i);
        [rxSamples,Parameters] = Channel(symbols_pulseShaped,Parameters);%+ (1/sqrt(2)) * (randn(size(symbols_pulseShaped)) + 1i * randn(size(symbols_pulseShaped))) * sigma(i);
        %rxSamples((sps*numTaps)/2+1+offset:sps:(length(symbols))*sps + (sps*numTaps)/2 +offset)
        % matched-filtering
        if(Parameters.PerfectChannelEst == "YS")
            % Matched-filtering
            rxSamplesMF = conv(rxSamples .* (conj(Parameters.Channel)./(abs(Parameters.Channel))),Parameters.h);
            rxSamplesMFOpt = rxSamplesMF(Parameters.delay:Parameters.sps:(Ns-1)*Parameters.sps + Parameters.delay);
        else
            % Matched-filtering
            rxSamplesMF = conv(rxSamples,Parameters.h);
            rxSamplesMFOpt = rxSamplesMF(Parameters.delay:Parameters.sps:(Ns-1)*Parameters.sps + Parameters.delay);
            % FFT Channel Estimation
            rxSamplesPilots = conj(Parameters.pilots.') .* rxSamplesMFOpt(Parameters.pilotVector);
            pilotChfft = fft(rxSamplesPilots,Np);
            pilotCh    = [pilotChfft(1:(Np-1)/2 + 1) zeros(1, Ns-Np) pilotChfft((Np-1)/2 + 2 : Np)];
            chEst      = ifft(pilotCh,Ns) * Ns/Np;
            rxSamplesMFOpt =  rxSamplesMFOpt .* conj(chEst)./(abs(chEst));
            %

        end
        rxSamplesEq =  rxSamplesMFDelayCorr;
        rxBits = MyDetectPSK(rxSamplesEq.',M);
        Nerrs(i) = Nerrs(i) + biterr(txBits,rxBits);
    end

end

BER_sim = Nerrs/(numBlocks*N);
semilogy(EbNodB,BER_sim,'LineWidth',1.5,'Marker','diamond');
hold on;
ber_theory = berfading(EbNodB,"psk",M,1);
%ber_theory = berawgn(EbNodB,"QAM",M);
semilogy(EbNodB,ber_theory,'LineWidth',1.5,'LineStyle',':');
grid on;
legend('simulated', 'theoretical')
xticks(EbNodB(1):2:EbNodB(end))
xlim([0 15])
ylim([10^(-5) 1])
xlabel('EbNodB')
ylabel('BER')
title('BPSK Modulation Scheme Simulation, # of symbols = 2')


% txBits = randi([0 1],N,1);
% [txBits,symbols] = MyPSK(txBits,M);
% %pulse shaping
% [symbols_pulseShaped, h, Eg] = PulseShape(symbols, sps, numTaps, f0, rolloff, pulseShape);


% bw = 1;
% figure;
% symbols_PS = symbols_pulseShaped(1:1024*sps);
% symbols_PS_F = fftshift(fft(symbols_PS));
% plot((-4096:4096-1)*bw/8192,10*log10(abs(symbols_PS_F)/max(abs(symbols_PS_F))),'LineWidth',1.25);
% xlim([-0.4 0.4])
% hold on;

% sps = 8;
% Ts = sps;
% alpha = 0.25;
% B0 = 1/Ts;
% f1 = (1-alpha)*B0;
%
% f = - 3:0.001:3;
%
% H = zeros(length(f),1);
%
% for i = 1:1:length(f)
%
% %     if(f(i)>-f1 && f(i)<f1)
% %         H(i) = sqrt(1/(2*B0));
% %     elseif(f(i) > -2*B0 + f1 && f(i) < 2*B0 - f1 )
% %         H(i) = sqrt(0.5 * (1/(2*B0)) * (1 + cos(pi/2 * ((abs(f(i))-f1)/(B0-f1)))));
% %     else
% %         H(i) = 0;
% %     end
%
%     if(f(i)>-f1 && f(i)<f1)
%         H(i) = sqrt(1/(2*B0));
%     else
%         H(i) = sqrt(0.5 * (1/(2*B0)) * (1 + cos(pi/2 * ((abs(f(i))-f1)/(B0-f1)))));
%     end
% end




% plot(f,10*log10(H/max(H)),'LineWidth',1.25,'LineStyle','--')
% legend('Simulated','Theoretical')
% grid on
% xlim([-0.4 0.4])
% xlabel('f in hz');
% ylabel('H(f) in log scale')
% title('Spectrum Plot')
% ylim([-30 5])