clear all
clc
USING_MATLAB = 0;

if (USING_MATLAB == 0)
  pkg load statistics
  pkg load communications
  pkg load signal
end

binary_data = [0, 1];
% No of sub carriers
N = 256;
% Desired data rate in bits/sec
Rd = 10^(6);
% Data rate of each sub carrier
Rs = Rd/N;
%symbol duration in sec calc
Tu = 2/Rs ;
% Freq separation
delF = 1/Tu;
% Number of timeslots to observe the simulation
Nt = 1;
% Random input data generation

P = zeros(1,199);
iter = zeros(1,199);
errors = zeros(1,199);
SNR = zeros(1,199);
sigma = 1;
ind = 0;
% Energy of each symbol 
for Es = 0.1:5:50 
%     ind = int16((Es-0.1)/0.05 + 1);
    ind = ind + 1;
    while(iter(ind)<1)
        iter(ind) = iter(ind)+1;
        if (USING_MATLAB) 
          data1 = randsample([0 1],N*Nt*2,true);  
          data2 = randsample([0 1],N*Nt*2,true);          
        else 
          rand_indices = randi(2, 1, N*Nt*2); 
          data1 = binary_data(rand_indices);
          rand_indices = randi(2, 1, N*Nt*2); 
          data2 = binary_data(rand_indices);
        end

        sigma = 0.1;
        % QPSK symbols - input to IFFT
        X1 = zeros(1,N*Nt);
        X2 = zeros(1,N*Nt);
        % QPSK symbol generation
        count = 0;
        for i=1:length(data1)/2
            count = count+1;
            switch bi2de([data1(2*i-1) data1(2*i)])
                case 0
                    X1(count) = complex(sqrt(Es/2),sqrt(Es/2));
                case 1
                    X1(count) = complex(-sqrt(Es/2),sqrt(Es/2));
                case 2
                    X1(count) = complex(-sqrt(Es/2),-sqrt(Es/2));
                case 3
                    X1(count) = complex(sqrt(Es/2),-sqrt(Es/2));
            end
        end
        count = 0;
        for i=1:length(data2)/2
            count = count+1;
            switch bi2de([data2(2*i-1) data2(2*i)])
                case 0
                    X2(count) = complex(sqrt(Es/2),sqrt(Es/2));
                case 1
                    X2(count) = complex(-sqrt(Es/2),sqrt(Es/2));
                case 2
                    X2(count) = complex(-sqrt(Es/2),-sqrt(Es/2));
                case 3
                    X2(count) = complex(sqrt(Es/2),-sqrt(Es/2));
            end
        end
        % Inverse FFT
        x1n = ifft(X1);
        x2n = ifft(X2);
        %n=0:1023;
        % 3-tap system
        taps=4;
        h = complex(random('norm',0,1,2,taps),random('norm',0,1,2,taps));
        % adding cyclic prefix
        for i=1:taps-1
            temp1(i) = x1n(length(x1n)-i+1);
            temp2(i) = x2n(length(x2n)-i+1);
        end
        x1n = [temp1 x1n];
        x2n = [temp2 x2n];
        % Channel response when xn is transmitted
        noise = random('norm',0,sigma,2,N);
        %noise = zeros(2,N);
        temp3 = conv(h(1,:),x1n) + conv(h(2,:),x2n);
        temp4 = conv(h(1,:),-fliplr(conj(x2n))) + conv(h(2,:),fliplr(conj(x1n)));
        y1n = temp3(taps:length(temp3)-taps+1) + noise(1,:);
        y2n = temp4(taps:length(temp4)-taps+1) + noise(2,:);
        % alamouti decoding begins :(
        H1 = fft(h(1,:),N);
        H2 = fft(h(2,:),N);
        Y1 = fft(y1n);
        Y2 = fft(y2n);
        % Estimate calculation
        U1 = (conj(H1).*Y1 + H2.*conj(Y2))./(abs(H1).^2 + abs(H2).^2);
        U2 = conj(H2).*Y1 - H1.*conj(Y2)./(abs(H1).^2 + abs(H2).^2);
        data_decode1=zeros(1,2*N);
        data_decode2=zeros(1,2*N);
    
        for i=1:N
            if(angle(U1(i))>=0 && angle(U1(i))< pi/2)
                data_decode1(2*i-1)=0;
                data_decode1(2*i)=0;
            elseif(angle(U1(i))>=pi/2 && angle(U1(i))< pi)
                data_decode1(2*i-1)=1;
                data_decode1(2*i)=0;
            elseif(angle(U1(i))>=-pi && angle(U1(i))< -pi/2)
                data_decode1(2*i-1)=0;
                data_decode1(2*i)=1;
            else
                data_decode1(2*i-1)=1;
                data_decode1(2*i)=1;
            end
        end
        for i=1:N
            if(angle(U2(i))>=0 && angle(U2(i))< pi/2)
                data_decode2(2*i-1)=0;
                data_decode2(2*i)=0;
            elseif(angle(U2(i))>=pi/2 && angle(U2(i))< pi)
                data_decode2(2*i-1)=1;
                data_decode2(2*i)=0;
            elseif(angle(U2(i))>=-pi && angle(U2(i))< -pi/2)
                data_decode2(2*i-1)=0;
                data_decode2(2*i)=1;
            else
                data_decode2(2*i-1)=1;
                data_decode2(2*i)=1;
            end
        end
        errors(ind) = errors(ind) + nnz(data_decode1-data1) + nnz(data_decode2-data2);  
        if(errors(ind)>=1200)
            break;
        end
    end
    P(ind) = errors(ind)/(2*N*iter(ind));
    SNR(ind) = 10*log10(Es/(sigma^2));
end

semilogy(SNR, P, "s", "linewidth", 2, "linestyle","-");
grid on
xlabel('Signal to Noise ratio (in db)');
ylabel('Bit error rate in logarithmic scale');
title("Simulating the Alamouti scheme performance");