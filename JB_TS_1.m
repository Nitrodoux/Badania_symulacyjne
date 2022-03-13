%%Projekt 1
%% Wyznaczenie ci¹g³ego modelu transmitancyjnegoo
%Zapisanie cz³onu opóŸniaj¹cego
[l,m]=pade(0.8,3)
%% Wyznaczenie modelu zastêpczego
%Zapis Gx(s)
Gx1=tf(l,m)
Gx2=tf([1],[2 1])
Gx=series(Gx1,Gx2)
%Zapis Gy(s)
Gy=tf([4],[1 0.5])
%Zapis Gz(s)
Gz=tf([3],[2 0 1])
G1=parallel(Gz,Gx)
G2=series(Gy,Gx)
% 1. Transmitancja wypadkowa
G=feedback(G1,G2)
%% Wyznaczenie minimalnej reprezentacji modelu zastêpczego
Gr=minreal(G)
%% Zapisanie minimalnej reprezentacji modelu w przestrzeni stanów
[liczl,mianl]=tfdata(Gr,'v')
%Macierze uk³adu w postaci jawnej
[Ar,Br,Cr,Dr]=tf2ss(liczl,mianl)
%Model obiektu w przestrzeni stanu
Gr_ss=ss(Ar,Br,Cr,Dr)
%% Wyznaczenie dyskretnego modelu metod¹ eksploratora zerowego rzêdu
Gd=c2d(Gr,0.1,'zoh')
%% Wyznaczenie równania charakterystycznego
I=eye(size(Ar)); %utworzenie macierzy jednostkowej I
s=sym('s');
Row_ch=det(s*I-Ar)% wyznacznik 
%% Wyznaczenie wartoœci w³asnych (biegunów)i zer uk³adu
bieguny=pole(Gr)
zera=zero(Gr)
%[V,d]=eig(Ar) to samo co pole
%% Mapa biegunów i zer
figure(1)
pzmap(Gr)
title('Mapa biegunów i zer obiektu');
xlabel('Oœ rzeczywista');
ylabel('Oœ urojona');
%% Czêstotliwoœci w³asne uk³adu, bezwymiarowy wspó³czynnik t³umienia
[wn z p ]=damp(Gr)
%czesto_wn – wektor czêstotliwoœci w³asnych uk³adu 
%wsp_tlum – wektor wspó³czynników t³umienia uk³adu
%p – wektor wartoœci w³asnych
%% OdpowiedŸ obiektu na wymuszenie impulsowe i skokowe jednostkowe
figure(2)
impulse(Gr,50)
xlabel('Czas [s]')
ylabel('Amplituda [m]')
title('OdpowiedŸ obiektu na wymuszenie impulsowe')

figure(3)
step(Gr,50,'k')
xlabel('Czas [s]')
ylabel('Amplituda [m]')
title('OdpowiedŸ obiektu na wymuszenie skokowe')
legend('Matlab')
%% Charakterystyka Nyquista i Bodego
figure(4)
str=nyquistoptions;
str.ShowFullContour='off';%wy³¹czenie ujemnych czêstotliwoœci
str.XLabel.String='Oœ rzeczywista';
str.YLabel.String='Oœ urojona';
str.Title.FontSize=12;
str.FreqUnits='Hz';%zmiana jednostek na osi x
nyquist(Gr,str);
title('Charakterystka Nyquista');

figure(5)
str1=bodeoptions;
str1.XLabel.String='Oœ rzeczywista';
str1.Ylabel.String='Oœ urojona';
str1.Ylabel.String={'Modu³' 'Faza'};
str1.Xlabel.String='Czêstotliwoœc';
str1.Title.FontSize=12;
str1.FreqUnits='Hz';
bode(Gr,str1);
title('Charakterystyka Bodego');
grid on;
%% OdpowiedŸ na wymuszenie sinusoidalne
t=[0:0.1:20]
fz=0.5*2*pi %czêstotliwoœæ f=0.5 Hz
az=1 %amplituda=1m
z1=az*sin(fz*t)
figure(6)
lsim(Gr,z1,t)
xlabel('Czas [s]')
ylabel('Amplituda [m]')
title('OdpowiedŸ obiektu na wymuszenie sinusoidalne')
%% Zapas fazy i modu³u
[Gm,Pm,Wgm,Wpm] = margin(Gr)
%Gm- wartoœc zapasu modu³u, Wgm- czêstotliwoœæ zapasu modu³u
%Pm- wartoœc zapasu fazy, Wpm- czêstotliwoœæ zapasu fazy
%% Linie pierwiastkowe i próba badañ stabilnoœci metod¹ Evansa
figure(7)
h = rlocusplot(Gr);
pp = getoptions(h); 
pp.Title.String = 'Linie pierwiastkowe';
pp.XLabel.String='Oœ rzeczywista';
pp.YLabel.String='Oœ urojona';
setoptions(h,pp); %Zatwierdzenie zmian na wykresie
%% Zapis uk³adu w przestrzeni stanu
Ar,Br,Cr,Dr
%% Sterowalnoœæ
S=ctrb(Ar,Br)
rank(Ar)
rank(S)
if rank(Ar)==rank(S)
    disp('Uk³ad jest sterowalny')
else
    disp('Uk³ad nie jest sterowalny')
end
%% Obserwowalnoœæ
O=obsv(Ar,Cr)
rank(Ar)
rank(O)
if rank(Ar)==rank(O)
    disp('Uk³ad jest obserwowalny')
else
    disp('Uk³ad nie jest obserwowalny')
end
%% Postacie kanonoiczne
%% Modalna/diagonalna
M=canon(Gr,'modal')
%% Obserwowalna=Sterowalna'
[Aster,Bster,Cster,Tster,kster]=ctrbf(Ar,Br,Cr)
[Aobsv,Bobsv,Cobsv,Tobsv,k] = obsvf(Ar,Br,Cr)
Ao=Aster'
Bo=Cster'
Co=Bster'
%canon(Gr,'companion')%postaæ sterowalna
%%
%[G,T]=minreal()
%Akalm=T*Ar*T'
%Bkalm=T*Br
%Ckalm=Cr*T
Gss=ss(Gr)
[sys,Ts]=canon(Gss,'modal')
Akalm=inv(Ts)*Ar*Ts
Bkalm=inv(Ts)*Br
Ckalm=Cr*Ts
%% Jordana
disp('Postaæ kanoniczna Jordana');
[V,D]=eig(Ar);
[V1,J1] = jordan(Ar);
jordan(Ar)
%% Wymuszenie skokowe simulnik
figure(8)
plot(sim.Time,sim.Data)
xlabel('Czas [s]');
ylabel('Amplituda [m]');
title('OdpowiedŸ na wymuszenie skokowe w simulinku')
%% K krytyczne
1.8445
for i=1.84:0.0001:1.85
    figure(9);
    hold on;
    step(feedback(i*Gr,1));
    hold on;
end

figure(10)
plot(kgran.Time,kgran.Data)
xlabel('Czas [s]');
ylabel('Amplituda [m]');
title('OdpowiedŸ uk³adu na wzmocnienie graniczne')
%}

