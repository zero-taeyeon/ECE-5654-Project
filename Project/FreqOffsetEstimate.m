function [FreqOffset]  =  FreqOffsetEstimate( InputSignal, MaxOffset, NumOffsets, h)

d = 2 * MaxOffset/(NumOffsets-1);
v = -MaxOffset:d:MaxOffset;
mfOutput = zeros(NumOffsets,length(InputSignal) + length(h)-1);
for i = 1:1:NumOffsets
    phase = exp(-1j*2*pi*v(i) * (0:1:length(InputSignal)-1));
    mfInput = phase .* InputSignal;
    mfOutput(i,:) = conv(mfInput,h);
end
% sps = 4;
% decisionFreq = sum((abs(mfOutput(:,length(PulseShape):sps:(100-1)*sps + length(PulseShape) ))).^2,2);
decisionFreq = sum((abs(mfOutput)).^2,2);
[~,freqEstId] = max(decisionFreq);
FreqOffset = v(freqEstId);

end