%%% Trabalho Final Vibrações%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

%%%%% AQUISIÇÃO DE DADOS SHAKER%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%Velocidade%%%%%%%%%%%%%%%%%%%%%%%%

VxT = xlsread('vel_GrupoA.xlsx');
V = VxT(:,2);
T = VxT(:,1);
for i = 1:700 
    V_(i) = V(i);
    T_(i) = T(i);
end


%%%%%%%% FORÇA%%%%%%%%%%%%%%%%%%%%%%%%

FxT = xlsread('forca_GrupoA.xlsx');
F = FxT(:,2);
Tf = FxT(:,1);
for i = 1:700 
    F_(i) = F(i);
    Tf_(i) = Tf(i);
end


%%%%%%% GRAFICO SHAKER%%%%%%%%%%%%%%%%%%%%%

%%%%%%%Calculo de c%%%%%%%%%%%%%%%%%%%%%%%%%

hold on
figure(1)
plot(F_,V_,'b')
title('Força-Velocidade')
grid on
axis([-120 80 -0.3 0.2])

hold off
[P3 P4] = ginput (2);
ceq = (P3(2)-P3(1))/(P4(2)-P4(1));
c = ceq/2;

%%%%%%%%%%%%%% AQUISIÇÃO DE DADOS IMPULSO%%%%%%%%%%%%%%%%%%%

resp(1).frontal=load('impulsofrontal01.txt')

hold on
figure(2)
plot(resp(1).frontal(:,1),resp(1).frontal(:,2),'g')
title('Resposta Livre')
axis([-0.5 2.5 -8*10^(-3) 8*10^(-3)])
hold off

TIME = resp(1).frontal(:,1);
U = resp(1).frontal(:,2);
M = 2; 
P = ginput(M);
d = (1/(M-1))*log(P(1,2)/P(M,2)); %Decremento logarítimo
% Coeficiente de amortecimento adimensional reduzido
Ksi = d/(2*pi);
% Determinando o período de oscilação
T=(P(M,1)-P(1,1))/(M-1); 
% Calculando a frequência natural amortecida
Wd = 2*pi/T;
% Calculando a frequência natural
Wn = Wd/sqrt(1-(Ksi^2));
% Equação para Senóide Amortecida
A = max(U)-min(U); %Amplitude do Sistema
x = A.*exp(-Ksi*Wn.*TIME).*cos(Wd.*TIME);

hold on
figure
plot(TIME,x)
title('Senóide Amortecida')
axis([-0.5 2.5 -6*10^(-3) 14*10^(-3)])
axis on
hold off

n = 2; 
E = ginput(n); 
p = polyfit(E(:,1),log(E(:,2)),1); 
Periodo = (E(n,1)-E(1,1))/(n-1);
Wn_n = sqrt((2*pi/Periodo)^2+(p(1,1))^2); %Frequencia natural
Ksi_n = abs(p(1,1)/Wn_n);


[nlin,ncol]=size(resp(1).frontal);
T=resp(1).frontal(2,1)-resp(1).frontal(1,1);

Fs = 1/T;            % Sampling frequency
L = nlin;             % Length of signal
t = (0:L-1)*T;        % Time vector

    X=resp(1).frontal(:,2);
    resp(1).Y = fft(X);
    resp(1).P2 = abs(resp(1).Y/L);
    resp(1).P1 = resp(1).P2(1:L/2+1);
    resp(1).P1(2:end-1) = 2*resp(1).P1(2:end-1);
    resp(1).f = Fs*(0:(L/2))/L;
    
hold on
figure(6)
plot(resp(1).f,resp(1).P1,'r')
title('FRF')
xlabel('f (Hz)')
ylabel('|Magnitude(m/s2)|')
axis([-5 120 -1*10^(-4) 6*10^(-4)])
grid on
grid minor
hold off


[w z]= ginput(3);

periodo2 = w(2) - w(1);

wd2 = w(3);
kissi2 = (w(2)-w(1))/wd2;

wn2 = wd2/sqrt(1-(kissi2^2));
MASSA = 272.5;

Keq = MASSA*wn2^2;
Kt = 120000;
Ks = Keq/(2 - Keq/Kt);

















