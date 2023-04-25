clear
clc
bits = [1 0 1 0 0 1 0 1 1 0];
qamConst = [-7 -5 -3 -1 1 3 5 7];
dataInMatrix = reshape(bits, 5, []); % columns are symbols
decdata = bit2int(dataInMatrix, 5);
symbol = [];
for i = 1:size(dataInMatrix,2)
    if decdata(i) == 0
        symbol(end+1) = -1-1j;
    elseif decdata(i) == 1
        symbol(end+1) = -3-1j;
    elseif decdata(i) == 2
        symbol(end+1) = -1-3j;
    elseif decdata(i) == 3
        symbol(end+1) = -3-3j; 

    elseif decdata(i) == 4
        symbol(end+1) = -1+1j;
    elseif decdata(i) == 5
        symbol(end+1) = -3+1j;
    elseif decdata(i) == 6
        symbol(end+1) = -1+3j;
    elseif decdata(i) == 7
        symbol(end+1) = -3+3j; 

    elseif decdata(i) == 8
        symbol(end+1) = 1-1j;
    elseif decdata(i) == 9
        symbol(end+1) = 3-1j;
    elseif decdata(i) == 10
        symbol(end+1) = 1-3j;
    elseif decdata(i) == 11
        symbol(end+1) = 3-3j; 

    elseif decdata(i) == 12
        symbol(end+1) = 1+1j;
    elseif decdata(i) == 13
        symbol(end+1) = 3+1j;
    elseif decdata(i) == 14
        symbol(end+1) = 1+3j;
    elseif decdata(i) == 15
        symbol(end+1) = 3+3j; 

    elseif decdata(i) == 16
        symbol(end+1) = -3-5j;
    elseif decdata(i) == 17
        symbol(end+1) = -5-1j;
    elseif decdata(i) == 18
        symbol(end+1) = -1-5j;
    elseif decdata(i) == 19
        symbol(end+1) = -5-3j; 

    elseif decdata(i) == 20
        symbol(end+1) = -3+5j;
    elseif decdata(i) == 21
        symbol(end+1) = -5+1j;
    elseif decdata(i) == 22
        symbol(end+1) = -1+5j;
    elseif decdata(i) == 23
        symbol(end+1) = -5+3j; 

    elseif decdata(i) == 24
        symbol(end+1) = 3-5j;
    elseif decdata(i) == 25
        symbol(end+1) = 5-1j;
    elseif decdata(i) == 26
        symbol(end+1) = 1-5j;
    elseif decdata(i) == 27
        symbol(end+1) = 5-3j;
     
    elseif decdata(i) == 28
        symbol(end+1) = 3+5j;
    elseif decdata(i) == 29
        symbol(end+1) = 5+1j;
    elseif decdata(i) == 30
        symbol(end+1) = 1+5j;
    elseif decdata(i) == 31
        symbol(end+1) = 5+3j;
    end
end

