%%Projekt 1
%% Wyznaczenie ci�g�ego modelu transmitancyjnegoo
%Zapisanie cz�onu op�niaj�cego
[l,m]=pade(0.8,3)
%% Wyznaczenie modelu zast�pczego
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
%% Wyznaczenie minimalnej reprezentacji modelu zast�pczego
Gr=minreal(G)
%% Zapisanie minimalnej reprezentacji modelu w przestrzeni stan�w
[liczl,mianl]=tfdata(Gr,'v')
%Macierze uk�adu w postaci jawnej
[Ar,Br,Cr,Dr]=tf2ss(liczl,mianl)
%Model obiektu w przestrzeni stanu
Gr_ss=ss(Ar,Br,Cr,Dr)
%% Wyznaczenie dyskretnego modelu metod� eksploratora zerowego rz�du
Gd=c2d(Gr,0.1,'zoh')
%% Wyznaczenie r�wnania charakterystycznego
I=eye(size(Ar)); %utworzenie macierzy jednostkowej I
s=sym('s');
Row_ch=det(s*I-Ar)% wyznacznik 
%% Wyznaczenie warto�ci w�asnych (biegun�w)i zer uk�adu
bieguny=pole(Gr)
zera=zero(Gr)
%[V,d]=eig(Ar) to samo co pole
%% Mapa biegun�w i zer
figure(1)
pzmap(Gr)
title('Mapa biegun�w i zer obiektu');
xlabel('O� rzeczywista');
ylabel('O� urojona');
%% Cz�stotliwo�ci w�asne uk�adu, bezwymiarowy wsp�czynnik t�umienia
[wn z p ]=damp(Gr)
%czesto_wn � wektor cz�stotliwo�ci w�asnych uk�adu 
%wsp_tlum � wektor wsp�czynnik�w t�umienia uk�adu
%p � wektor warto�ci w�asnych
%% Odpowied� obiektu na wymuszenie impulsowe i skokowe jednostkowe
figure(2)
impulse(Gr,50)
xlabel('Czas [s]')
ylabel('Amplituda [m]')
title('Odpowied� obiektu na wymuszenie impulsowe')

figure(3)
step(Gr,50,'k')
xlabel('Czas [s]')
ylabel('Amplituda [m]')
title('Odpowied� obiektu na wymuszenie skokowe')
legend('Matlab')
%% Charakterystyka Nyquista i Bodego
figure(4)
str=nyquistoptions;
str.ShowFullContour='off';%wy��czenie ujemnych cz�stotliwo�ci
str.XLabel.String='O� rzeczywista';
str.YLabel.String='O� urojona';
str.Title.FontSize=12;
str.FreqUnits='Hz';%zmiana jednostek na osi x
nyquist(Gr,str);
title('Charakterystka Nyquista');

figure(5)
str1=bodeoptions;
str1.XLabel.String='O� rzeczywista';
str1.Ylabel.String='O� urojona';
str1.Ylabel.String={'Modu�' 'Faza'};
str1.Xlabel.String='Cz�stotliwo�c';
str1.Title.FontSize=12;
str1.FreqUnits='Hz';
bode(Gr,str1);
title('Charakterystyka Bodego');
grid on;
%% Odpowied� na wymuszenie sinusoidalne
t=[0:0.1:20]
fz=0.5*2*pi %cz�stotliwo�� f=0.5 Hz
az=1 %amplituda=1m
z1=az*sin(fz*t)
figure(6)
lsim(Gr,z1,t)
xlabel('Czas [s]')
ylabel('Amplituda [m]')
title('Odpowied� obiektu na wymuszenie sinusoidalne')
%% Zapas fazy i modu�u
[Gm,Pm,Wgm,Wpm] = margin(Gr)
%Gm- warto�c zapasu modu�u, Wgm- cz�stotliwo�� zapasu modu�u
%Pm- warto�c zapasu fazy, Wpm- cz�stotliwo�� zapasu fazy
%% Linie pierwiastkowe i pr�ba bada� stabilno�ci metod� Evansa
figure(7)
h = rlocusplot(Gr);
pp = getoptions(h); 
pp.Title.String = 'Linie pierwiastkowe';
pp.XLabel.String='O� rzeczywista';
pp.YLabel.String='O� urojona';
setoptions(h,pp); %Zatwierdzenie zmian na wykresie
%% Zapis uk�adu w przestrzeni stanu
Ar,Br,Cr,Dr
%% Sterowalno��
S=ctrb(Ar,Br)
rank(Ar)
rank(S)
if rank(Ar)==rank(S)
    disp('Uk�ad jest sterowalny')
else
    disp('Uk�ad nie jest sterowalny')
end
%% Obserwowalno��
O=obsv(Ar,Cr)
rank(Ar)
rank(O)
if rank(Ar)==rank(O)
    disp('Uk�ad jest obserwowalny')
else
    disp('Uk�ad nie jest obserwowalny')
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
%canon(Gr,'companion')%posta� sterowalna
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
disp('Posta� kanoniczna Jordana');
[V,D]=eig(Ar);
[V1,J1] = jordan(Ar);
jordan(Ar)
%% Wymuszenie skokowe simulnik
figure(8)
plot(sim.Time,sim.Data)
xlabel('Czas [s]');
ylabel('Amplituda [m]');
title('Odpowied� na wymuszenie skokowe w simulinku')
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
title('Odpowied� uk�adu na wzmocnienie graniczne')
%}

