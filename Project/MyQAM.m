function [bits, symbols] = MyQAM(bits,numSymbols)

bitsPerSymbol = log2(numSymbols);

if(mod(length(bits),bitsPerSymbol) ~= 0)
    bits = [bits; zeros(bitsPerSymbol - mod(length(bits),bitsPerSymbol),1)];
end

switch bitsPerSymbol

    case 4        % 16QAM Modulation

        %         symbols = (1/sqrt(10)) * ( (1-2*bits(1:4:end)) .* (2 - (1-2*bits(3:4:end))) ...
        %             + 1i*(1-2*bits(2:4:end)) .* (2 - (1-2*bits(4:4:end))));

        symbols = (1/sqrt(10)) * ( (2*bits(1:4:end)-1) .* ((3-2*bits(2:4:end))) ...
            + 1i*(1-2*bits(3:4:end)) .* (3-2*bits(4:4:end)));

        % symbols = (1/sqrt(10)) * ( (1-2*bits(1:4:end)) .* (3 - 2*bits(2:4:end))) ...
        %              + 1i*(1-2*bits(3:4:end)) .* (3 - 2*bits(4:4:end)));

    case 6        % 64QAM Modulation

        %         re_part =  (1-2*bits(1:6:end)) .* (4 - (1-2*bits(3:6:end)) ...
        %             .* (2 - (1-2*bits(5:6:end))));

        re_part =  (2*bits(1:6:end)-1) .* (4 - (2*bits(2:6:end)-1) ...
            .* (3 -2*bits(3:6:end)));

        im_part =  (1-2*bits(4:6:end)) .* (4 - (2*bits(5:6:end)-1) ...
            .* (3 -2*bits(6:6:end)));

        symbols = (1/sqrt(42)) * (re_part + 1i*im_part);

    case 5        % 32QAM Modulation
        dataInMatrix = reshape(bits, 5, []); % columns are symbols
        decdata = bit2int(dataInMatrix, 5);
        symbols = zeros(length(decdata),1);
        %         for i = 1:size(dataInMatrix,2)
        %             if decdata(i) == 0
        %                 symbols(i) = -1-1j;
        %             elseif decdata(i) == 1
        %                 symbols(i) = -3-1j;
        %             elseif decdata(i) == 2
        %                 symbols(i) = -1-3j;
        %             elseif decdata(i) == 3
        %                 symbols(i) = -3-3j;
        %
        %             elseif decdata(i) == 4
        %                 symbols(i) = -1+1j;
        %             elseif decdata(i) == 5
        %                 symbols(i) = -3+1j;
        %             elseif decdata(i) == 6
        %                 symbols(i) = -1+3j;
        %             elseif decdata(i) == 7
        %                 symbols(i) = -3+3j;
        %
        %             elseif decdata(i) == 8
        %                 symbols(i) = 1-1j;
        %             elseif decdata(i) == 9
        %                 symbols(i) = 3-1j;
        %             elseif decdata(i) == 10
        %                 symbols(i) = 1-3j;
        %             elseif decdata(i) == 11
        %                 symbols(i) = 3-3j;
        %
        %             elseif decdata(i) == 12
        %                 symbols(i) = 1+1j;
        %             elseif decdata(i) == 13
        %                 symbols(i) = 3+1j;
        %             elseif decdata(i) == 14
        %                 symbols(i) = 1+3j;
        %             elseif decdata(i) == 15
        %                 symbols(i) = 3+3j;
        %
        %             elseif decdata(i) == 16
        %                 symbols(i) = -3-5j;
        %             elseif decdata(i) == 17
        %                 symbols(i) = -5-1j;
        %             elseif decdata(i) == 18
        %                 symbols(i) = -1-5j;
        %             elseif decdata(i) == 19
        %                 symbols(i) = -5-3j;
        %
        %             elseif decdata(i) == 20
        %                 symbols(i) = -3+5j;
        %             elseif decdata(i) == 21
        %                 symbols(i) = -5+1j;
        %             elseif decdata(i) == 22
        %                 symbols(i) = -1+5j;
        %             elseif decdata(i) == 23
        %                 symbols(i) = -5+3j;
        %
        %             elseif decdata(i) == 24
        %                 symbols(i) = 3-5j;
        %             elseif decdata(i) == 25
        %                 symbols(i) = 5-1j;
        %             elseif decdata(i) == 26
        %                 symbols(i) = 1-5j;
        %             elseif decdata(i) == 27
        %                 symbols(i) = 5-3j;
        %
        %             elseif decdata(i) == 28
        %                 symbols(i) = 3+5j;
        %             elseif decdata(i) == 29
        %                 symbols(i) = 5+1j;
        %             elseif decdata(i) == 30
        %                 symbols(i) = 1+5j;
        %             elseif decdata(i) == 31
        %                 symbols(i) = 5+3j;
        %             end
        %         end

        for i = 1:size(dataInMatrix,2)
            if decdata(i) == 0
                symbols(i) = -3+5j;
            elseif decdata(i) == 1
                symbols(i) = -1+5j;
            elseif decdata(i) == 2
                symbols(i) = -3-5j;
            elseif decdata(i) == 3
                symbols(i) = -1-5j;

            elseif decdata(i) == 4
                symbols(i) = -5+3j;
            elseif decdata(i) == 5
                symbols(i) = -5+1j;
            elseif decdata(i) == 6
                symbols(i) = -5-3j;
            elseif decdata(i) == 7
                symbols(i) = -5-1j;

            elseif decdata(i) == 8
                symbols(i) = -1+3j;
            elseif decdata(i) == 9
                symbols(i) = -1+1j;
            elseif decdata(i) == 10
                symbols(i) = -1-3j;
            elseif decdata(i) == 11
                symbols(i) = -1-1j;

            elseif decdata(i) == 12
                symbols(i) = -3+3j;
            elseif decdata(i) == 13
                symbols(i) = -3+1j;
            elseif decdata(i) == 14
                symbols(i) = -3-3j;
            elseif decdata(i) == 15
                symbols(i) = -3-1j;

            elseif decdata(i) == 16
                symbols(i) = +3+5j;
            elseif decdata(i) == 17
                symbols(i) = +1+5j;
            elseif decdata(i) == 18
                symbols(i) = 3-5j;
            elseif decdata(i) == 19
                symbols(i) = 1-5j;

            elseif decdata(i) == 20
                symbols(i) = 5+3j;
            elseif decdata(i) == 21
                symbols(i) = 5+1j;
            elseif decdata(i) == 22
                symbols(i) = 5-3j;
            elseif decdata(i) == 23
                symbols(i) = 5-1j;

            elseif decdata(i) == 24
                symbols(i) = 1+3j;
            elseif decdata(i) == 25
                symbols(i) = 1+1j;
            elseif decdata(i) == 26
                symbols(i) = 1-3j;
            elseif decdata(i) == 27
                symbols(i) = 1-1j;

            elseif decdata(i) == 28
                symbols(i) = 3+3j;
            elseif decdata(i) == 29
                symbols(i) = 3+1j;
            elseif decdata(i) == 30
                symbols(i) = 3-3j;
            elseif decdata(i) == 31
                symbols(i) = 3-1j;
            end
        end
        symbols = (1/sqrt(22)) * symbols;

    case 8  % 256QAM Modulation

        re_part = (1-2*bits(1:8:end)) .* (8 - (1-2*bits(3:8:end)) .* ...
            (4 - (1-2*bits(5:8:end)) .* (2 - (1-2*bits(7:8:end)))));

        im_part = (1-2*bits(2:8:end)) .* (8 - (1-2*bits(4:8:end)) .* ...
            (4 - (1-2*bits(6:8:end)) .* (2 - (1-2*bits(8:8:end)))));

        symbols = (1/sqrt(170)) * (re_part + 1i*im_part);
    otherwise
        disp('Invalid Modulation Order selected');


end


end