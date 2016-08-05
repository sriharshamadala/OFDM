clear all
clc
N=256
G1 = 23;% octal 7 corresponds to binary 111 n1 = m1 + m0 + m-1
G2 = 35;% octal 3 corresponds to binary 011 n1 = m0  + m-1
constLen = 5;   % Constraint length
% adding a tail of zeros for CC decoding
datainf = randsample([0 1],N-constLen,true);
datainf = [datainf zeros(1,constLen)];
% Create the trellis that represents the convolutional code
convCodeTrellis = poly2trellis(constLen, [G1 G2]);
codedWord = convenc(datainf, convCodeTrellis);
noise = random('norm',0,.1,1,2*N);
ncodedWord = codedWord + noise;

dd = quantiz(ncodedWord,[0 0.5 0.8],[3 2 1 0]);
data_decode = vitdec(dd,convCodeTrellis,5*constLen,'term','soft',2);

errors = nnz(data_decode - datainf)