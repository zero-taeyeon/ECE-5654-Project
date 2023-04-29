% This code evaluates the design of your transmitter/receiver
clear; clc; close all
%%% Parameters that can be changed by the student  %%%%
Parameters.ChannelType = 'AWGN'; % 'AWGN', 'RAYL', or 'RICE'
Parameters.fd = 100;              % max Doppler (10 or 100) only used for fading
Parameters.K = 10;               % Ricean K-factor; only used in Ricean fading
Parameters.PayloadSize = 100;     % this is completely up to the student
% This is the number of info bits per packet.
% The total packet size depends on the
% overhead included and coding rate.

Parameters.PerfectChannelEst = 'NO';% 'YS' - puts channel info into
% Parameters.Channel
% 'NO' - channel estimation must be done
%  separately

%%% Parameters that SHOULD NOT bbe changed by the student  %%%%
Parameters.NumSinusoids = 100;      % used by Channel function
Parameters.fs = 10e7;               % used by Channel function

% what I added

% Pulse shaping parameters
Parameters.sps = 4;                 % samples for symbol
Parameters.numTaps = 10;            % Filter taps
Parameters.rolloff   = 0.25;        % roll-off factor


Parameters.pulseShape = 'SqRa';     % pulse shape

% Channel Encoding
Parameters.encoder = 'Polar';
Parameters.numiter = 4;

%%% Frequency Synchronization
Parameters.freqOffset = 0.0;                                 % frequency offset
Parameters.MaxOffset = 0.15;               % maximum frequency offset
Parameters.NumOffsets= 1;                        % number of offsets

%%% timing offset params
maxOffset = 3;
Parameters.timeOffset = randi([1,maxOffset]);

%%% frame synchronization

% rangen of SNRs to test over
if Parameters.ChannelType == "AWGN"
    SNRdB = 0:8;
else
    SNRdB = 0:2:20;
end

% Number of packets to simulate to estimate performance.  You can change
% this for testing, but for final plots you should set at 100000
NumPackets = 10;
bit_errors = zeros(length(SNRdB),1);
correct_bits = zeros(length(SNRdB),1);

% Modulation Parameters
ModType = ["PSK","PSK","QAM","QAM","QAM"];
M = [2,4,16,32,64];


for j = 1:1:5
    Parameters.Modtype = ModType(j);         % modulation type
    Parameters.M = M(j);                  % modulation level

    for i=1:length(SNRdB)
        tic
        Parameters.SNR = 10.^(SNRdB(i)/10);
        bit_errors(i) = 0;
        correct_bits(i) = 0;

        for k=1:NumPackets

            b = round(rand(1,Parameters.PayloadSize));

            [OutputSamples,Parameters] = MyTransmitter(b, Parameters);

            [ReceivedSamples, Parameters] = Channel(OutputSamples, Parameters);

            %NOTE: ReceivedSamples should have 1 row per Rx antenna and one
            % column per sample.

            b_est = MyReceiver(ReceivedSamples, Parameters).';

            bit_errors(i) = bit_errors(i) + sum(abs(b-b_est));
            if sum(abs(b-b_est)) == 0
                correct_bits(i) = correct_bits(i) + Parameters.PayloadSize;
            end

        end

        BER(i,j) = bit_errors(i)/Parameters.PayloadSize/NumPackets;

        Throughput(i,j) = correct_bits(i)/(size(ReceivedSamples,2)*(1/Parameters.fs));

        % Display the elapsed time
        if((i) == 1)
            toc;
            elapsed_time = toc;
            fprintf('Time taken to complete one SNR: %.2f seconds, total time: %.2f seconds\n', elapsed_time,elapsed_time*length(SNRdB));
        end
    end

end
Legend = {'BPSK','QPSK','16QAM','32QAM','64QAM'};
% BER plot
figure
for i = 1:1:5
    semilogy(SNRdB, BER(:,i))
    hold on;
end
legend(Legend)
grid on;
xlabel('SNR (dB)')
ylabel('BER')
xlim([SNRdB(1) SNRdB(end)])
ylim([10^(-4) 1])

% Througput Plot
figure
for i = 1:1:5
    plot(SNRdB, Throughput(:,i)/100)
    hold on;
end
legend(Legend)
grid on;
xlabel('SNR (dB)')
ylabel('Througput (b/s)')
xlim([SNRdB(1) SNRdB(end)])
%ylim([10^(-4) 1])

save("AWGN_perfect.mat", 'BER', 'Throughput')
