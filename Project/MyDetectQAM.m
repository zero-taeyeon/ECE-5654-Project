function rxbits = MyDetectQAM(samples,M)

bitsPerSymbol = log2(M);
ip = 0:M-1;
op = bitxor(ip,floor(ip/2));
d = dictionary(ip+1,op);

b_gry = zeros(M,bitsPerSymbol);
for i = bitsPerSymbol:-1:1
    b_gry(:,i) = mod(op,2);
    op = floor(op/2);
end
refbits = reshape(b_gry.',size(b_gry,1)*bitsPerSymbol,1);
[~,refSymbols] = MyQAM(refbits,M);


[~ , idx] = min(abs(samples - refSymbols.'),[],2);
gry_dec = d(idx);
b_gry = zeros(length(samples),bitsPerSymbol);
for i = bitsPerSymbol:-1:1
    b_gry(:,i) = mod(gry_dec,2);
    gry_dec = floor(gry_dec/2);
end
rxbits = reshape(b_gry.',size(b_gry,1)*bitsPerSymbol,1);

% switch M
%     case 16
% 
%         %         refSymbols = 1/sqrt(10) *  [ 1 + 1i  1 + 3i  3 + 3i  3 + 1i  ...
%         %             3 - 1i  3 - 3i  1 - 3i  1 - 1i ...
%         %             -1 - 1i -1 - 3i -3 - 3i -3 - 1i ...
%         %             -3 + 1i -3 + 3i -1 + 3i -1 + 1i];
% 
%         [~ , idx] = min(abs(samples - refSymbols.'),[],2);
%         gry_dec = d(idx);
%         b_gry = zeros(length(samples),bitsPerSymbol);
%         for i = bitsPerSymbol:-1:1
%             b_gry(:,i) = mod(gry_dec,2);
%             gry_dec = floor(gry_dec/2);
%         end
%         rxbits = reshape(b_gry.',size(b_gry,1)*bitsPerSymbol,1);
% 
%     case 64
% 
%         %         refSymbols = 1/sqrt(10) *  [ 1 + 1i  1 + 3i  3 + 3i  3 + 1i  ...
%         %             3 - 1i  3 - 3i  1 - 3i  1 - 1i ...
%         %             -1 - 1i -1 - 3i -3 - 3i -3 - 1i ...
%         %             -3 + 1i -3 + 3i -1 + 3i -1 + 1i];
% 
% 
%     otherwise



end