% Channel code
function [OutputSamples, Parameters] = Channel(InputSamples, Parameters) 

% noise standard deviation - assumes that four samples are used per symbol,
% but that the noise bandwidth is equal to the symbol rate.
sigma = sqrt(4/2/Parameters.SNR);

normalized_samples = InputSamples/sqrt(mean(abs(InputSamples.^2)));

if Parameters.ChannelType ~= 'AWGN'
    f = Parameters.fd*sin(2*pi*rand(1,Parameters.NumSinusoids));
    Theta = 2*pi*rand(1,Parameters.NumSinusoids);
    t = Parameters.fs*(0:(length(InputSamples)-1));
    tmp = 0;
    for ii=1:length(f)
        tmp = tmp + exp(1i*(2*pi*f(ii)*t+Theta(ii)));
    end
    Parameters.Channel = tmp/sqrt(Parameters.NumSinusoids);
else
    Parameters.Channel = ones(1,length(InputSamples));
end

if Parameters.ChannelType == 'RICE'
    Parameters.Channel = sqrt(Parameters.K/(Parameters.K+1))*exp(1i*2*pi*rand) + ...
                         sqrt(1/(Parameters.K+1))*Parameters.Channel;
end


OutputSamples = Parameters.Channel.*InputSamples + ...
                sigma*(randn(1,length(InputSamples))+...
                1i*randn(1,length(InputSamples)));


