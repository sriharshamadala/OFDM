clear all
clc

% Only needed for octave
pkg load statistics
pkg load communications
pkg load signal

number_of_subcarriers = 256;

% Desired data rate in bits/sec
Rd = 10^(6);

% Data rate of each sub carrier
Rs = Rd/number_of_subcarriers;

% Symbol duration in sec calc
Tu = 2/Rs ;

% Freq separation
delF = 1/Tu;

% Number of timeslots to observe the simulation
Nt = 1;

% Random input data generation
sigma = 1;
SNR = [0:4:32 35:3:50];
Es = (sigma^2).*10.^(SNR./10);
P = zeros(1,length(Es));
iter = zeros(1,length(Es));
errors = zeros(1,length(Es));
binary_data = [0, 1];

% Energy of each symbol 
for ind=1:length(Es) 
    while(iter(ind)<5000)
        iter(ind) = iter(ind)+1;
        %data = randsample([0 1],number_of_subcarriers*Nt*2,true);          
        rand_indices = randi(2, 1, number_of_subcarriers*Nt*2); 
        data = binary_data(rand_indices);
        
        % QPSK symbols - input to IFFT
        X = zeros(1,number_of_subcarriers*Nt);
        
        % QPSK symbol generation
        count = 0;
        for i=1:length(data)/2
            count = count+1;
            switch bi2de([data(2*i-1) data(2*i)])
                case 0
                    X(count) = complex(sqrt(Es(ind)/2),sqrt(Es(ind)/2));
                case 1
                    X(count) = complex(-sqrt(Es(ind)/2),sqrt(Es(ind)/2));
                case 2
                    X(count) = complex(-sqrt(Es(ind)/2),-sqrt(Es(ind)/2));
                case 3
                    X(count) = complex(sqrt(Es(ind)/2),-sqrt(Es(ind)/2));
            end
        end
        
        % 4-tap system
        taps=4;
        h = complex(random('norm',0,1,1,taps),random('norm',0,1,1,taps));
        H = fft(h,number_of_subcarriers);
        
        % Precoding
        %X_precod = X./abs(H);

        % Inverse FFT
        xn = ifft(X);         
        
        % adding cyclic prefix
        %for i=1:taps-1
        %     temp1(i) = xn(length(xn)-i+1);
        %end
        %xn = [temp1 xn];

        % Channel response when xn is transmitted
        noise = (1/sqrt(2))*complex(random('norm',0,sigma,1,number_of_subcarriers),...
                                      random('norm',0,sigma,1,number_of_subcarriers));
        %temp2 = cconv(h,xn,number_of_subcarriers);
        
        % Uncomment the previous line in matlab as circular convolution is not implemented in octave.
        temp2 = ifft(H.*X);
        
        %yn = temp2(taps:length(temp2)-taps+1);% + noise;     
        yn = temp2 + noise;
        Y = fft(yn);
        Y = Y./H;
        %Y = Y.*complex(cos(angle(H)),-sin(angle(H)));
        data_decode=zeros(1,2*number_of_subcarriers);
        for i=1:number_of_subcarriers
            if(angle(Y(i))>=0 && angle(Y(i))< pi/2)
                data_decode(2*i-1)=0;
                data_decode(2*i)=0;
            elseif(angle(Y(i))>=pi/2 && angle(Y(i))< pi)
                data_decode(2*i-1)=1;
                data_decode(2*i)=0;
            elseif(angle(Y(i))>=-pi && angle(Y(i))< -pi/2)
                data_decode(2*i-1)=0;
                data_decode(2*i)=1;
            else
                data_decode(2*i-1)=1;
                data_decode(2*i)=1;
            end
        end
        errors(ind) = errors(ind) + nnz(data_decode-data);  
        if(errors(ind)>=1000)
            break;
        end
    end
    P(ind) = errors(ind)/(2*number_of_subcarriers*iter(ind));
end

semilogy(SNR, P, "s", "linewidth", 2, "linestyle","-");
xlabel('Signal to Noise ration (in db)');
ylabel('Bit error rate in logarithmic scale');
grid on
title(sprintf('BER vs SNR of OFDM in a %d-tap system with gaussian noise', taps))