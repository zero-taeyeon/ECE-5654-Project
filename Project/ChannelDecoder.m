% Channel Encoder
function [decodedBits, Parameters] = ChannelDecoder(rxBitsllr, Parameters)

if(Parameters.encoder == "Polar")
    turbodec = comm.TurboDecoder(Parameters.trellis,Parameters.intrlvrIndices,Parameters.numiter);
    decodedBits = turbodec(-rxBitsllr.');
elseif(Parameters.encoder == "LDPC")
% 
% 
%     k = K(r);                           % number of bits of message per block
%     R = rate(r);                        % coded rate
%     H = dvbs2ldpc(R);
%     cfgLDPCEnc = ldpcEncoderConfig(H);
%     cfgLDPCDec = ldpcDecoderConfig(H);
end


end