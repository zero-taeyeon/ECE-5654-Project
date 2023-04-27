function [bits, symbols] = MyPSK(bits,numSymbols)


bitsPerSymbol = log2(numSymbols);

if(mod(length(bits),bitsPerSymbol) ~= 0)
    bits = [bits; zeros(bitsPerSymbol - mod(length(bits),bitsPerSymbol),1)];
end

%% PSK Encoding

b_gry = reshape(bits,bitsPerSymbol, length(bits)/bitsPerSymbol);
b_gry = b_gry.'; % N x log2*(M) matrix
b_bin = zeros(size(b_gry));
b_bin(:,1) = b_gry(:,1);

for i = 2:bitsPerSymbol
    b_bin(:,i) = bitxor(b_bin(:,i-1),b_gry(:,i)); % matrix after gray coding
end

dec_mat = 2.^(bitsPerSymbol-1:-1:0);

b_dec = sum(dec_mat .* b_bin,2); % decimal vector

symbols = exp(1i * ((2*pi/numSymbols) * b_dec ));

end