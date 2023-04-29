function [y,h,Eg] = PulseShape( x, Ns, N, rolloff, PulseShape)

if(1)
sps = Ns;
% Ts = sps;            
% t = -floor(N*sps/2):ceil(N*sps/2);

switch PulseShape
    case 'SQAR'
        h = ones(1,sps);
    case 'SINC'
        h = sinc(-numTaps/2:1/Ns:numTaps/2);
    case 'RaCo'
        %h = (1/Ts) .* sinc(t/Ts) .* (cos(pi * rolloff * t/Ts)./(1 - (2*rolloff*t/Ts).^2));
        h = rcosdesign(rolloff,N, Ns);
    case 'SqRa'
        h = rcosdesign(rolloff,N, Ns, 'sqrt');
    otherwise
        disp('Pulse shape not supported')
end

Eg = sum(h.^2);
h = (1/sqrt(Eg))*h;
X = zeros(length(x),sps);
X(:,1) = x;
X = reshape(X.',1,length(x)*sps);
y =  conv(X,h);

else
    y = x.';
    h = 1;
    Eg = 1;
end

end



