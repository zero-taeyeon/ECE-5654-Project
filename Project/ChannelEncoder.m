% Channel Encoder
function [OutputBits, Parameters] = ChannelEncoder(InputBits, Parameters)

if(Parameters.encoder == "Polar")

    K = Parameters.PayloadSize;
    trellis = poly2trellis(4,[13 15 17],13);
    n = log2(trellis.numOutputSymbols);
    numTails = log2(trellis.numStates)*n;
    Parameters.E = K*(2*n - 1) + 2*numTails;           % Output codeword packet length
    Parameters.rate = K/Parameters.E; % Coding rate
    intrlvrIndices = randperm(K);
    turboenc = comm.TurboEncoder(trellis,intrlvrIndices);
    Parameters.trellis = trellis;
    Parameters.intrlvrIndices = intrlvrIndices;
    %turbodec = comm.TurboDecoder(trellis,intrlvrIndices,numiter);
    OutputBits = turboenc(InputBits.');
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