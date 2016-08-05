clear all
clc
USING_MATLAB = 0;

nf_subcarriers = 256;
min_nf_errors = 100;
max_iters = 1000;
binary_data = [0, 1];

% Desired data rate in bits/sec
Rd = 10^(6);

% Data rate of each sub carrier
Rs = Rd/nf_subcarriers;

% Symbol duration in sec calc
Tu = 2/Rs ;

% Number of timeslots to observe the simulation
Nt = 1;

% Random input data generation
SNR = [0:4:32 35:3:50];
sigma = 1;
Es = (sigma^2).*10.^(SNR./10);
P = zeros(1,length(Es));
iter = zeros(1,length(Es));
errors = zeros(1,length(Es));

% Energy of each symbol 
for ind = 1:length(Es) 
  while (iter(ind) <= max_iters)
    iter(ind) = iter(ind) + 1;           
    G1 = 23;        % octal 7 corresponds to binary 111 n1 = m1 + m0 + m-1 
    G2 = 35;        % octal 3 corresponds to binary 011 n1 = m0  + m-1 
    constLen = 5;   % Constraint length  
    % adding a tail of zeros for CC decoding
    if (USING_MATLAB)
      datainf = randsample([0 1], nf_subcarriers-constLen, true);
    else
      rand_indices = randi(2, 1, nf_subcarriers-constLen);
      datainf = binary_data(rand_indices);
    end
    datainf = [datainf zeros(1, constLen)]; 
    
    % Create the trellis that represents the convolutional code
    convCodeTrellis = poly2trellis(constLen, [G1 G2]);
    codedWord1 = convenc(datainf, convCodeTrellis);
    data = codedWord1;
            
    % QPSK symbols - input to IFFT
    X = zeros(1, nf_subcarriers*Nt);
    % QPSK symbol generation
    count = 0;
    for i = 1:length(data)/2
      count = count+1;
      switch bi2de([data(2*i-1) data(2*i)])
        case 0
          X(count) = complex(-sqrt(Es(ind)/2), -sqrt(Es(ind)/2));
        case 1
          X(count) = complex(sqrt(Es(ind)/2), -sqrt(Es(ind)/2));
        case 2
          X(count) = complex(-sqrt(Es(ind)/2), sqrt(Es(ind)/2));
        case 3
          X(count) = complex(sqrt(Es(ind)/2), sqrt(Es(ind)/2));
      end
    end
    
    % Interleaving block
    s2 = length(X);
    interShift = s2/64;
    matrix = reshape(X, interShift, 64);
    interleaveData = matrix';
    Xint = interleaveData(:)';
    
    % Inverse FFT
    xn = ifft(Xint);
    
    %n=0:1023;
    taps=4;
    h = complex(random('norm',0,1,1,taps),random('norm',0,1,1,taps));
    H = fft(h, nf_subcarriers);

%         % adding cyclic prefix
%         for i=1:taps-1
%             temp1(i) = xn(length(xn)-i+1);
%         end
%         xn = [temp1 xn];
    % Channel response when xn is transmitted
    noise = (1/sqrt(2))*complex(random('norm',0,sigma,1,nf_subcarriers),...
                                 random('norm',0,sigma,1,nf_subcarriers));
    if (USING_MATLAB)                             
      temp2 = cconv(h,xn,nf_subcarriers);
    else
      temp2 = ifft(H.*Xint);
    end

    yn = temp2 + noise;        
    Y = fft(yn);
    Y = Y./H;
   
    %Deinterleaving starts 
    data_mod = reshape(Y,64,interShift);       
    deintlvddata=data_mod';
    deintlvddata=deintlvddata(:)';
    
    % Viterbi decoder
    ddR = quantiz(real(deintlvddata/sqrt(Es(ind))),[-0.4 0 0.4],[3 2 1 0]);
    ddI = quantiz(imag(deintlvddata/sqrt(Es(ind))),[-0.4 0 0.4],[3 2 1 0]);
    dd = zeros(1,2*length(ddR));
    for i = 1:length(ddR)
        dd(2*i-1) = ddR(i);
        dd(2*i) = ddI(i);
    end
    receivedBits = vitdec(dd,convCodeTrellis,5*constLen,'term','soft',2);
%         receivedBits = vitdec(data_decode,convCodeTrellis,5*constLen,'term','hard');
%         Y = complex(receivedBitsReal,receivedBitsImg);
%         Y = receivedBits;
    data_decode=receivedBits; 
%         for i=1:length(receivedBits)
%             data_decode(2*i-1) = receivedBitsReal(i);
%             data_decode(2*i) = receivedBitsImg(i);
%         end
    
    errors(ind) = errors(ind) + nnz(data_decode-datainf);   
    if(errors(ind)>=min_nf_errors)
        break;
    end
  end
  P(ind) = errors(ind)/(nf_subcarriers*iter(ind));
end