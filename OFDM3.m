clear all
clc
USING_MATLAB = 0;

if (USING_MATLAB == 0)
  pkg load statistics
  pkg load communications
  pkg load signal
end

nf_subcarriers = 256;
min_nf_errors = 1000;
max_iters = 5000;

% Desired data rate in bits/sec
Rd = 10^(6);

% Data rate of each sub carrier
Rs = Rd/nf_subcarriers;

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
for ind = 1:length(Es) 
  while (iter(ind) < max_iters)
    iter(ind) = iter(ind) + 1;
    if (USING_MATLAB) 
      data = randsample([0 1], nf_subcarriers*Nt*2, true);          
    else 
      rand_indices = randi(2, 1, nf_subcarriers*Nt*2); 
      data = binary_data(rand_indices);
    end
    
    % QPSK symbols - input to IFFT
    X = zeros(1, nf_subcarriers*Nt);
    
    % QPSK symbol generation
    count = 0;
    for i=1:length(data)/2
      count = count+1;
      switch bi2de([data(2*i-1) data(2*i)])
        case 0
          X(count) = complex(sqrt(Es(ind)/2), sqrt(Es(ind)/2));
        case 1
          X(count) = complex(-sqrt(Es(ind)/2), sqrt(Es(ind)/2));
        case 2
          X(count) = complex(-sqrt(Es(ind)/2), -sqrt(Es(ind)/2));
        case 3
          X(count) = complex(sqrt(Es(ind)/2), -sqrt(Es(ind)/2));
      end
    end
    
    % 4-tap system
    taps=4;
    h = complex(random('norm', 0, 1, 1, taps),random('norm', 0, 1, 1, taps));
    H = fft(h, nf_subcarriers);
    
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
    noise = (1/sqrt(2))*complex(random('norm', 0, sigma, 1, nf_subcarriers),...
                                  random('norm', 0, sigma, 1, nf_subcarriers));
    
    if (USING_MATLAB)    
      temp2 = cconv(h, xn, nf_subcarriers);
    else 
      % Circular convolution is not implemented in octave. Using covolution theorem for DFT.
      temp2 = ifft(H.*X);
    end
    
    %yn = temp2(taps:length(temp2)-taps+1);% + noise;     
    yn = temp2 + noise;
    Y = fft(yn);
    Y = Y./H;
    %Y = Y.*complex(cos(angle(H)),-sin(angle(H)));
    data_decode = zeros(1, 2*nf_subcarriers);
    for i=1:nf_subcarriers
      if ((angle(Y(i)) >= 0) && (angle(Y(i)) < pi/2))
        data_decode(2*i-1) = 0;
        data_decode(2*i) = 0;
      elseif ((angle(Y(i)) >= pi/2) && (angle(Y(i)) < pi))
        data_decode(2*i-1) = 1;
        data_decode(2*i) = 0;
      elseif ((angle(Y(i)) >= -pi) && (angle(Y(i))< -pi/2))
        data_decode(2*i-1) = 0;
        data_decode(2*i) = 1;
      else
        data_decode(2*i-1) = 1;
        data_decode(2*i) = 1;
      end
    end
    errors(ind) = errors(ind) + nnz(data_decode-data);  
    if (errors(ind) >= min_nf_errors)
      % We need atleast min_nf_errors errors for each SNR to rely on the probability of error computation.
      break;
    end
  end
  P(ind) = errors(ind)/(2*nf_subcarriers*iter(ind));
end

semilogy(SNR, P, "s", "linewidth", 2, "linestyle","-");
grid on
xlabel('Signal to Noise ratio (in db)');
ylabel('Bit error rate in logarithmic scale');
title(sprintf('BER vs SNR of OFDM in a %d-tap system with gaussian noise', taps));