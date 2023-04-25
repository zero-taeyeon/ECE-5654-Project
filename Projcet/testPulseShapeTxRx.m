%% Start of the code
clc;clear;close all;

%% carrier
chan.f                      = 6e6;                              % Carrier frequency
chan.lambda                 = 3e8/chan.f;                       % Carrier wavelength
%% Antenna Structure and Spacing
chan.Nt                     = 1;                                % Transmit Antennae
chan.Mr                     = 1;                                % Receive Antennae
chan.d_lambda_ratio         = 10;                               % (d/lambda) ratio
chan.rxantennaSpacing       = chan.d_lambda_ratio*chan.lambda;  % Rx Antenna Spacing
chan.txantennaSpacing       = chan.d_lambda_ratio*chan.lambda;  % Tx Antenna Spacing
chan.rxarray                = [0:chan.rxantennaSpacing:(chan.Mr -1)*chan.rxantennaSpacing];
chan.txarray                = [0:chan.txantennaSpacing:(chan.Nt -1)*chan.txantennaSpacing];
%% Multipath, Rayleigh and Ricean Effects, Doppler
chan.riceanK_factor    = 0;                                     % Ricean Factor
chan.Np                = 1;                                     % Discrete path delays % Flat Fading for Np = 1
chan.fs                = 10e6;                                  % Sampling frequency
chan.Ts                = 1/chan.fs;                             % Sampling Time
chan.N                 = 25;                                   % Un-resolvable multipath components
chan.c                 = 1/sqrt(chan.N);                        % Scaling Factor
chan.fd                = 100;                                   % Maximum Doppler Frequency (Hz)
%% Time to Simulate
chan.time       = 1;                                           % Time to Simulate
chan.Ns         = chan.fs * chan.time;                          % Number of samples to simulate
chan.t          = [0:chan.Ts:chan.time-chan.Ts];                % Time vector
%% Angle Spread, Angle of Arrival & Angle of Departure
% Creates Angle's with Uniform Distribution with min value and max values
% as Input parameters
% For creating Gaussian distribution, use @normrnd iwith mean and variance
% in radians
% Ex: chan.theta1 = AngleDistribution(@normrand, mean , anglespread, [chan.N 1]);
% Angle of arrival of multipath component relative to mobile velocity
chan.theta1 = AngleDistribution(@unifrnd, 0,2*pi,[chan.N 1]);
% Initial Angle of arrival due to minor difference in distance travelled due to unresolvable multipath components
chan.theta2 = AngleDistribution(@unifrnd, 0,2*pi,[chan.N 1]);
% Angle of departure at the Tx Antenna
chan.theta3 = AngleDistribution(@unifrnd, 0,2*pi,[chan.N 1]);    


%% pulse shape parameters
sps = 4;            % samples for symbol
numTaps = 25;
rolloff   = 0.25;
pulseShape = 'RaCo';
% Assuming sample period is 0.1 micro sec, so sample rate is 10 Msps, 
% so *symbol* period is sps

if(pulseShape == "SQAR")
    delay = sps;
else
    delay = sps*numTaps+1;%length(h);
end

%% modulation Params
% generate bits
M = 2;
EbNodB = 0:15 ;
R = 1; % uncoded BPSK (1 bit/symbol)
SNR = EbNodB + 10*log10(log2(M));
EbNo = 10.^(EbNodB/10);
sigma = sqrt((10.^(-SNR/10)));
N = chan.fs/sps; % number of bits of message per block
numBlocks = 5;
Nerrs = zeros(size(EbNodB));
offset = 0;

for j = 1:1:numBlocks

    txBits = randi([0 1],N,1);
    % encoding
    
    % bits to symbols
    [txBits,symbols] = MyPSK(txBits,M);
    % pulse shaping
    [symbols_pulseShaped, h, Eg] = PulseShape(symbols, sps, numTaps, rolloff, pulseShape);

    for i = 1:1:length(EbNo)
        % add fading
        % mimo_H = Channel(chan);
        % add noise
        rxSamples = symbols_pulseShaped + (1/sqrt(2)) * (randn(size(symbols_pulseShaped)) + 1i * randn(size(symbols_pulseShaped))) * sigma(i);
        
        
        % frequency synchronizatiom
        % time synchrnonization
        % phase synchronization
        rxSamplesMF = conv(rxSamples,h);
        rxSamplesMFDelayCorr = rxSamplesMF(delay+offset:sps:(N-1)*sps + delay +offset);
        rxBits = MyDetectPSK(rxSamplesMFDelayCorr.',M);
        Nerrs(i) = Nerrs(i) + biterr(txBits,rxBits);
    end

end

BER_sim = Nerrs/(numBlocks*N);
semilogy(EbNodB,BER_sim,'LineWidth',1.5,'Marker','diamond');
hold on;
ber_theory = berawgn(EbNodB,"psk",M,'nondiff');
semilogy(EbNodB,ber_theory,'LineWidth',1.5,'LineStyle',':');
grid on;
legend('simulated', 'theoretical')
xticks(EbNodB(1):2:EbNodB(end))
xlim([0 15])
ylim([10^(-5) 1])
xlabel('EbNodB')
ylabel('BER')
title('BPSK Modulation Scheme Simulation, # of symbols = 2')
