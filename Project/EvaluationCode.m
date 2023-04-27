% This code evaluates the design of your transmitter/receiver
clear all
%%% Parameters that can be changed by the student  %%%%
Parameters.ChannelType = 'AWGN'; % 'AWGN', 'RAYL', or 'RICE'
Parameters.fd = 10;              % max Doppler (10 or 100) only used for fading
Parameters.K = 10;               % Ricean K-factor; only used in Ricean fading
Parameters.PayloadSize = 100;    % this is completely up to the student
% This is the number of info bits per packet.
% The total packet size depends on the
% overhead included and coding rate.

Parameters.PerfectChannelEst = 'YS';% 'YS' - puts channel info into
% Parameters.Channel
% 'NO' - channel estimation must be done
% separately
% what I added
Parameters.sps = 4;            % samples for symbol
Parameters.numTaps = 25;
Parameters.rolloff   = 0.25;

Parameters.Modtype = 'PSK';         % modulation type
Parameters.M = 2;                   % modulation level
Parameters.pulseShape = 'SqRa';     % pulse shape
% Channel Encoding
Parameters.encoder = 'Polar';
Parameters.numiter = 4;

%%% Parameters that SHOULD NOT bbe changed by the student  %%%%
Parameters.NumSinusoids = 100;      % used by Channel function
Parameters.fs = 1e7;               % used by Channel function



%%% Frequency Synchronization
Parameters.freqOffset = 0.0;                                 % frequency offset
Parameters.MaxOffset = 0.15;               % maximum frequency offset
Parameters.NumOffsets= 200;                        % number of offsets


% rangen of SNRs to test over
if Parameters.ChannelType == "AWGN"
    SNRdB = 0:8;
else
    SNRdB = 0:9;
end

% Number of packets to simulate to estimate performance.  You can change
% this for testing, but for final plots you should set at 100000
NumPackets = 1000;
SNRdB = 6;

for i=1:length(SNRdB)

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

    BER(i) = bit_errors(i)/Parameters.PayloadSize/NumPackets;

    Throughput(i) = correct_bits(i)/(size(ReceivedSamples,2)*Parameters.fs);

end

% BER plot
figure
semilogy(SNRdB, BER)
xlabel('SNR (dB)')
ylabel('BER')

% Througput Plot
figure
plot(SNRdB, Throughput)
xlabel('SNR (dB)')
ylabel('Througput (b/s)')

