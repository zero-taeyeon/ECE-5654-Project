function rxbits = MyDetectPSK(samples,M)

refSymbols = exp(1i * ((2*pi/M) * (0:M-1) ) );

[~ , idx] = min(abs(samples - refSymbols).^2,[],2);

%rx_b_bin = dec2bin(idx-1);
bitsPerSymbol = log2(M);
b_bin = zeros(length(samples),bitsPerSymbol);
dec = idx - 1;
for i = bitsPerSymbol:-1:1
    b_bin(:,i) = mod(dec,2);
    dec = floor(dec/2);
end
b_gry = zeros(size(b_bin));
b_gry(:,1) = b_bin(:,1);
for i = 2:bitsPerSymbol
    b_gry(:,i) = bitxor(b_bin(:,i-1),b_bin(:,i));
end

rxbits = reshape(b_gry.',size(b_bin,1)*bitsPerSymbol,1);
end