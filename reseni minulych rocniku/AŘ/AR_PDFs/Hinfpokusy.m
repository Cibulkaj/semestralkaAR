%Ukazkove priklady navrhu pres mixed sensitivity problem

%System 1.radu
%-------------
%Navrhove pozadavky - sirka pasma 10 rad/s
%tlumeni NF a VF poruch

s=zpk('s');
G=(5)/(s+1);
MS=1.05; %max.peak citlivostni funkce
AS=0.1; %utlum poruch na nizkych frekvencich
WS=20; %sirka pasma

W1=(s/MS+WS)/(s+AS*WS);
%bode(1/W1)
%nula vahove funkce urcuje sirku pasma regulace, pol zase miru potlaceni NF poruch, 
%nezavisi na zesileni W1 - to se kompenzuje hodnotou GAMMA, kterou algoritmus nalezne 
%(musi platit Bodeho integralni formule)
W2=1;
W3=1;

[K,CL,GAM]=mixsyn(G,W1,[],[]);
L=G*K; S=inv(1+L); T=1-S;
sigma(tf(S),'g',tf(T),'r',GAM/W1,'g-.')
%regulator ani system nema integracni slozku, vychazi lead-lag kompenzator
%s jednim rychlym polem jdoucim k - nekonecnu


%pridani penalizace rizeni konstantou W2
MS=1.05; %max.peak citlivostni funkce
AS=0.1; %utlum poruch na nizkych frekvencich
WS=20; %sirka pasma

W1=(s/MS+WS)/(s+AS*WS);
W2=0.00001  %velmi mala penalizace, omezeni neni aktivni a neovlivni tvar Tjw
%volba konstantni fce W2 omezuje sirku pasma, pri zachovani spodni hranice potlaceni NF poruch
[K,CL,GAM]=mixsyn(G,W1,W2,[]);
L=G*K; S=inv(1+L); T=1-S;
sigma(S,'g',T,'r',K/S,'m',GAM*G/ss(W2),'r-.',GAM/W1,'g-.')

%pridani penalizace rizeni konstantou W2
MS=1.05; %max.peak citlivostni funkce
AS=0.1; %utlum poruch na nizkych frekvencich
WS=20; %sirka pasma

W1=(s/MS+WS)/(s+AS*WS);
W2=1.5  %velka penalizace rizeni
%snizeni sirky pasma a horsi potlaceni poruch
[K,CL,GAM]=mixsyn(G,W1,W2,[]);
L=G*K; S=inv(1+L); T=1-S;
sigma(S,'g',T,'r',K/S,'m',GAM*G/ss(W2),'r-.',GAM/W1,'g-.')



%SOUCASNA PENALIZACE Sjw a Tjw - tvarovani Ljw vymezenim zakazanych oblasti
%VYJDE REGULATOR 3.RADU, KTERY LZE VYPUSTENIM RYCHLYCH POLU A NULY
%REDUKOVAT NA LEAD-LEG KOMPENZATOR 2.54*(S+1)/(S+0.5)
s=zpk('s');
G=(5)/(s+1); 
MS=1.05; %max.peak citlivostni funkce
AS=0.1; %utlum poruch na nizkych frekvencich
WS=5; %sirka pasma
W1=(s/MS+WS)/(s+AS*WS);
MT=1.05; %max peak kompl. citl. fce
AT=.005; %utlum na vysokych frekvencich
WT=50; %sirka pasma
W3=(s+WT/MT)/(AT*s+WT);
close all;
%bode(W1);
%hold;
%bode(1/W3);
[K0,CL,GAM]=mixsyn(G,W1,[],W3);
%K0 je optimalni regulator pro rozsireny system Gext
L=G*K0; S=inv(1+L); T=1-S;
figure;
sigma(S,'g',T,'r',GAM/W3,'r-.',GAM/W1,'g-.')
figure;
%Tvarovani Ljw
%W1 omezuje L zdola na nizkych frekvencich (odregulovani poruch a sledovani
%ref.signalu), 1/W3 pak shora na vysokych (robustnost)
sigma(L,W1,1/W3)


%TRETI RAD
%------------

%POZADAVEK NA PRIDANI INTEGRATORU
%PROBLEM - NULA REGULATORU VELMI BLIZKO K NULE KTERA KRATI ZAMYSLENY NULOVY
%POL 
s=zpk('s');
G=(5)/(s+1)^3; 
Gext=G/s; %rozsireni systemu o integrator
MS=1.05; %max.peak citlivostni funkce
AS=0.1; %utlum poruch na nizkych frekvencich
WS=5; %sirka pasma
W1=(s/MS+WS)/(s+AS*WS);
MT=1.05; %max peak kompl. citl. fce
AT=.05; %utlum na vysokych frekvencich
WT=20; %sirka pasma
W3=(s+WT/MT)/(AT*s+WT);
%bode(1/W3)
[K0,CL,GAM]=mixsyn(Gext,W1,[],W3);
%K0 je optimalni regulator pro rozsireny system Gext
L=Gext*K0; S=inv(1+L); T=1-S;
sigma(S,'g',T,'r',GAM/W3,'r-.',GAM/W1,'g-.')
%Resenim pro puvodni system G je K0/s 
K=tf(K0)/s;
L=G*K; S=inv(1+L); T=1-S;
sigma(S,'g',T,'r',GAM/W3,'r-.',GAM/W1,'g-.')

%PRIDANI INTEGRATORU - OPRAVA K1
s=zpk('s');
G=(5)/(s+1)^3; 
Gext=G/s; %rozsireni systemu o integrator
W1=1;
V=(s+10)/s; %presun polu otevrene smycky z nuly do -10 pozad.sirka pasma 10 rad/s
sigma(1/(1*V)); %chovani Sjw pro nizke frekvence
W2=0.0001; %mala penalizace rizeni 
P=ss([W1*V W1*Gext;0 W2;-V -Gext]); %vstupy porucha vystup, rizeni, vystupy W1*y,W2*u,vstup regulatoru
P=minreal(P);
P=mktito(P,1,1);
[K,CL,GAM,INFO] = hinfsyn(P);
%kompenzator pro puvodni problem
K0=K/s;
L=G*K0; S=inv(1+L); T=1-S;
close all;
sigma(S,'g',T,'r',GAM/V,'g-.',GAM*G/ss(W2),'r-.')
zpk(K0)
CL=L/(1+L)
zpk(CL)
pole(CL)
%VYSLEDNY REGULATOR MA INTEGRATOR, UZAVRENA SMYCKA OBSAHUJE ZAMYSLENY POL V
%-10, NENI TO OVSEM DOMINANTNI POL VLIVEM MALE PENALIZACE RIZENI

%PRIDANI PENALIZACE NA RIZENI - dosazena sirka pasma 1 rad/s
%kompenzator 3.rad + integrator
V=(s+10)/s;
W1=1;
W2=0.001;
P=ss([W1*V W1*Gext;0 W2;-V -Gext]);
P=minreal(P);
P=mktito(P,1,1);
[K,CL,GAM,INFO] = hinfsyn(P);
%kompenzator pro puvodni problem
K0=K/s;
L=G*K0; S=inv(1+L); T=1-S;
close all;
sigma(S,'g',T,'r',GAM/(W1*V),'g-.',GAM*G/(W2*V),'r-.')
zpk(K0)
CL=L/(1+L)
zpk(CL)
pole(CL)
%omezeni na citlivostni fci rizeni
U=K0/(1+L);
figure;
sigma(U,GAM/ss(V*W2))
%REGULATOR MA VF POKLES 20db/dek vlivem integratoru, pri pozadavku na vetsi
%nutno modifikovat fci W2 na c(1+rs) a prislusne upravit celou ulohu kvuli
%resitelnosti stavove realizace


%JESTE VETSI PENALIZACE
W2=10;
P=ss([W1*V W1*Gext;0 W2;-V -Gext]);
P=minreal(P);
P=mktito(P,1,1);
[K,CL,GAM,INFO] = hinfsyn(P);
%kompenzator pro puvodni problem
K0=K/s;
L=G*K0; S=inv(1+L); T=1-S;
close all;
sigma(S,'g',T,'r',GAM/V,'g-.',GAM*G/ss(W2),'r-.')
zpk(K0)
CL=L/(1+L)
zpk(CL)
pole(CL)
%omezeni na citlivostni fci rizeni
U=K0/(1+L);
figure;
sigma(U,GAM/ss(V*W2))


%PRESUNUTI DVOU POLU K2
s=zpk('s');
G=(5)/(s+1)^3; 
Gext=G/s; %rozsireni systemu o integrator
W1=1;
V=(s^2+2*10*s+100)/(s*(s+1)); %presun dvou polu otevrene smycky z nuly do -10 pozad.sirka pasma 10 rad/s
%PRI MALE PENALIZACI RIZENI KONSTANTOU W2 BUDOU DVA POLY V -10 DOMINANTNI
sigma(1/(1*V)); %chovani Sjw pro nizke frekvence
W2=0.01; %mala penalizace rizeni 
P=ss([W1*V W1*Gext;0 W2;-V -Gext]);
P=minreal(P);
P=mktito(P,1,1);
[K,CL,GAM,INFO] = hinfsyn(P);
%kompenzator pro puvodni problem
K0=K/s;
L=G*K0; S=inv(1+L); T=1-S;
close all;
sigma(S,'g',T,'r',GAM/V,'g-.',GAM*G/ss(W2),'r-.')
zpk(K0)
CL=L/(1+L)
zpk(CL)
pole(CL)

%PRESUNUTI ctyr POLU K3
s=zpk('s');
G=(5)/(s+1)^3; 
Gext=G/s; %rozsireni systemu o integrator
W1=1;
V=((s^2+2*10*s+100)*(s+50)^2)/(s*(s+1)^3); %presun ctyr polu otevrene smycky na dva dominantni do -10 a dva do -20 pozad.sirka pasma 10 rad/s

%Korekce pozadavku na citl.funkci filtrem W1
%W1=(s+1)^3/(s*(s+50)^2);
W1=1;

%PRI MALE PENALIZACI RIZENI KONSTANTOU W2 BUDOU DVA POLY V -10 DOMINANTNI
sigma(1/(V*W1)); %chovani Sjw pro nizke frekvence
W2=0.1; %mala penalizace rizeni 
P=ss([W1*V W1*Gext;0 W2;-V -Gext]);
P=minreal(P);
P=mktito(P,1,1);
[K,CL,GAM,INFO] = hinfsyn(P);
%kompenzator pro puvodni problem
K0=K/s;
L=G*K0; S=inv(1+L); T=1-S;
close all;
sigma(S,'g',T,'r',GAM/(V*W1),'g-.',GAM*G/ss(W2),'r-.')
zpk(K0)
CL=L/(1+L)
zpk(CL)
pole(CL)







%3.RAD SE SPECIFIKACI W1 A W3
%PROBLEM - REGULATOR NELZE ROZUMNYM ZPUSOBEM REDUKOVAT, MA VELMI RYCHLE
%MODY, NAVIC NEMA INTEGRATOR
s=zpk('s');
G=(5)/((s+1)*(s+2)*(s+3)); 
MS=1.05; %max.peak citlivostni funkce
AS=0.1; %utlum poruch na nizkych frekvencich
WS=5; %sirka pasma
W1=(s/MS+WS)/(s+AS*WS);
MT=1.05; %max peak kompl. citl. fce
AT=.005; %utlum na vysokych frekvencich
WT=50; %sirka pasma
W3=(s+WT/MT)/(AT*s+WT);
close all;
bode(W1);
hold;
bode(1/W3);
[K0,CL,GAM]=mixsyn(G,W1,[],W3);
%K0 je optimalni regulator pro rozsireny system Gext
L=G*K0; S=inv(1+L); T=1-S;
figure;
sigma(S,'g',T,'r',GAM/W3,'r-.',GAM/W1,'g-.')
figure;
sigma(L,W1,1/W3)





%NEMIN.FAZOVY SYSTEM S 1 NESTABILNI NULOU
%SPATNY NAVRH, KTERY NERESPEKTUJE NESTABILNU NULU
%PODLE BODEHO INTEGRALNICH OMEZENI NELZE DOSAHNOUT TAKTO VELKE SIRKY PASMA
s=zpk('s');
G=(5)*(s-1)/((s+1)*(s+2)); 
MS=1.05; %max.peak citlivostni funkce
AS=0.1; %utlum poruch na nizkych frekvencich
WS=20; %sirka pasma
W1=(s/MS+WS)/(s+AS*WS);
MT=1.05; %max peak kompl. citl. fce
AT=.005; %utlum na vysokych frekvencich
WT=50; %sirka pasma
W3=(s+WT/MT)/(AT*s+WT);
close all;
bode(W1);
hold;
bode(1/W3);
[K0,CL,GAM]=mixsyn(G,W1,[],W3);
%K0 je optimalni regulator pro rozsireny system Gext
L=G*K0; S=inv(1+L); T=1-S;
figure;
sigma(S,'g',T,'r',GAM/W3,'r-.',GAM/W1,'g-.')
figure;
sigma(L,W1,1/W3)












%NEMIN.FAZOVY SYSTEM S 1 NESTABILNI NULOU
%omezeni sirky pasma, ted uz spravne
%regulator ovsem nema integrator
s=zpk('s');
G=(5)*(s-1)/((s+1)*(s+2)); 
MS=1.05; %max.peak citlivostni funkce
AS=0.1; %utlum poruch na nizkych frekvencich
WS=0.1; %sirka pasma
W1=(s/MS+WS)/(s+AS*WS);
MT=1.05; %max peak kompl. citl. fce
AT=.005; %utlum na vysokych frekvencich
WT=1; %sirka pasma
W3=(s+WT/MT)/(AT*s+WT);
close all;
bode(W1);
hold;
bode(1/W3);
[K0,CL,GAM]=mixsyn(G,W1,[],W3);
%K0 je optimalni regulator pro rozsireny system Gext
L=G*K0; S=inv(1+L); T=1-S;
figure;
sigma(S,'g',T,'r',GAM/W3,'r-.',GAM/W1,'g-.')
figure;
sigma(L,W1,1/W3)












%NEMIN.FAZOVY SYSTEM S 1 NESTABILNI NULOU
%omezeni sirky pasma
%doplneni o integrator


s=zpk('s');
G=(5)*(-s+1)/((s+1)*(s+2)); 
Gext=G/s; %rozsireni systemu o integrator
W1=1;
V=(s^2+2*s+1)/(s*(s+1)); %presun nuloveho polu a polu -1 na 2x -1

%Korekce pozadavku na citl.funkci filtrem W1
%W1=(s+1)^3/(s*(s+50)^2);
W1=1;

%PRI MALE PENALIZACI RIZENI KONSTANTOU W2 BUDOU DVA POLY V -10 DOMINANTNI
sigma(1/(V*W1)); %chovani Sjw pro nizke frekvence
W2=0.1; %mala penalizace rizeni 
P=ss([W1*V W1*Gext;0 W2;-V -Gext]);
P=minreal(P);
P=mktito(P,1,1);
[K,CL,GAM,INFO] = hinfsyn(P);
%kompenzator pro puvodni problem
K0=K/s;
L=G*K0; S=inv(1+L); T=1-S;
close all;
sigma(S,'g',T,'r',GAM/(V*W1),'g-.',GAM*G/ss(W2),'r-.')
zpk(K0)
CL=L/(1+L)
zpk(CL)
pole(CL)
%Velky peak citlivostni funkce, je treba snizit sirku pasma


s=zpk('s');
G=(5)*(-s+1)/((s+1)*(s+2)); 
Gext=G/s; %rozsireni systemu o integrator
W1=1;
V=(s^2+s+0.25)/(s*(s+1)); %presun nuloveho polu a polu -1 na 2x -1
%Korekce pozadavku na citl.funkci filtrem W1
%W1=(s+1)^3/(s*(s+50)^2);
W1=1;
%PRI MALE PENALIZACI RIZENI KONSTANTOU W2 BUDOU DVA POLY V -10 DOMINANTNI
sigma(1/(V*W1)); %chovani Sjw pro nizke frekvence
W2=2; %mala penalizace rizeni 
P=ss([W1*V W1*Gext;0 W2;-V -Gext]);
P=minreal(P);
P=mktito(P,1,1);
[K,CL,GAM,INFO] = hinfsyn(P);
%kompenzator pro puvodni problem
K0=K/s;
L=G*K0; S=inv(1+L); T=1-S;
close all;
sigma(S,'g',T,'r',GAM/(V*W1),'g-.',GAM*G/ss(W2),'r-.')
CL=L/(1+L)
zpk(CL)
pole(CL)
zpk(K0)




%NESTABILNI SYSTEM - 1 nestabilni pol

%SPATNY NAVRH - REGULATOR KRATI NESTABILNI POL NESTABILNI NULOU
%Pri minimalni chybe v modelu systemu bude uzavr. smycka nestabilni !!!

s=zpk('s');
G=(5)*(s+1)/((s-5)*(s+2)); 
Gext=G/s; %rozsireni systemu o integrator
W1=1;
%V=(s^2+2*s+1)/(s*(s+2)); %presun nuloveho polu a polu -2 na 2x -1
V=(s+100)/s;
%Korekce pozadavku na citl.funkci filtrem W1
%W1=(s+1)^3/(s*(s+50)^2);
W1=1;
%PRI MALE PENALIZACI RIZENI KONSTANTOU W2 BUDOU DVA POLY V -10 DOMINANTNI
sigma(1/(V*W1)); %chovani Sjw pro nizke frekvence
W2=0.01; %mala penalizace rizeni 
P=ss([W1*V W1*Gext;0 W2;-V -Gext]);
P=minreal(P);
P=mktito(P,1,1);
[K,CL,GAM,INFO] = hinfsyn(P,'display','on');
%kompenzator pro puvodni problem
K0=K/s;
L=G*K0; S=inv(1+L); T=1-S;
close all;
sigma(S,'g',T,'r',GAM/(V*W1),'g-.',GAM*G/ss(W2),'r-.')
CL=L/(1+L)
zpk(CL)
pole(CL)
zpk(K0)



%POKUS C.2, PRESUN PROBLEMATICKEHO NESTAB.POLU JINAM
%ZNOVU SPATNY NAVRH VLIVEM BODEHO INTEGRALNICH OMEZENI - SIRKA PASMA MUSI
%BYT VETSI NEZ HODNOTA NESTABILNIHO POLU

s=zpk('s');
G=(5)*(s+1)/((s-5)*(s+2)); 
Gext=G/s; %rozsireni systemu o integrator
W1=1;
V=(s+0.1)^2/(s*(s+5)); %presun nuloveho polu a polu -2 na 2x -1
%Korekce pozadavku na citl.funkci filtrem W1
%W1=(s+1)^3/(s*(s+50)^2);
W1=1
W1=(s+5)/(s+0.2);
%PRI MALE PENALIZACI RIZENI KONSTANTOU W2 BUDOU DVA POLY V -10 DOMINANTNI
sigma(1/(V*W1)); %chovani Sjw pro nizke frekvence
W2=1; %mala penalizace rizeni 
sigma(1/(V*W2)); %chovani Ujw pro vysoke frekvence
sigma(G/(V*W2)); %chovani Tjw pro vusoke f
P=ss([W1*V W1*Gext;0 W2;-V -Gext]);
P=minreal(P);
P=mktito(P,1,1);
[K,CL,GAM,INFO] = hinfsyn(P,'display','on');
%kompenzator pro puvodni problem
K0=K/s;
L=G*K0; S=inv(1+L); T=1-S;
close all;
sigma(S,'g',T,'r',GAM/(V*W1),'g-.',GAM*G/ss(V*W2),'r-.')
CL=L/(1+L)
zpk(CL)
pole(CL)
zpk(K0)

%POKUS C.3 - SPRAVNY NAVRH, SIRKA PASMA JE VETSI NEZ HODNOTA NESTABILNIHO
%POLU
s=zpk('s');
G=(5)*(s+1)/((s-5)*(s+2)); 
Gext=G/s; %rozsireni systemu o integrator
W1=1;
V=(s+50)^2/(s*(s-5)); %presun nuloveho polu a polu -2 na 2x -1
%Korekce pozadavku na citl.funkci filtrem W1
%W1=(s+1)^3/(s*(s+50)^2);
W1=1
W1=(s-5)/(s+50);
%PRI MALE PENALIZACI RIZENI KONSTANTOU W2 BUDOU DVA POLY V -10 DOMINANTNI
sigma(1/(V*W1)); %chovani Sjw pro nizke frekvence
W2=0.005; %mala penalizace rizeni 
%W2=100*(s*(s-5))/(s+10000)^2;%modifikace W2 pro penalizaci rizeni na VF, vysledny prenos 1/(W2*V) je 1/(s+50)^2, tedy pozadavek na utlum 40db/dek od 50 rad/s
W2=(s+60)/(s+50000); %modifikace W2 pro VF roll-off
sigma(1/(V*W2)); %chovani Ujw pro vysoke frekvence
sigma(G/(V*W2)); %chovani Tjw pro vusoke f
P=ss([W1*V W1*Gext;0 W2;-V -Gext]);
P=minreal(P);
P=mktito(P,1,1);
[K,CL,GAM,INFO] = hinfsyn(P,'display','on');
%kompenzator pro puvodni problem
K0=K/s;
L=G*K0; S=inv(1+L); T=1-S;
close all;
sigma(S,'g',T,'r',GAM/(V*W1),'g-.',GAM*G/ss(V*W2),'r-.')
CL=L/(1+L)
zpk(CL)
pole(CL)
zpk(K0)
figure;
sigma(K0/(1+L),'r',GAM/(V*W2),'r-.');
%Vysledny regulator lze redukovat na PID regulator (2 nuly,1 pol v nule, druhy stabilni)
%nebo PID + LP filter pro VF roll-off






%Kmitavy system se slabe tlumenym modem, pozadavek na vyssi sirku pasma
s=zpk('s');
G=1/(s^2+2*0.1*s+1); 
Gext=G/s; %rozsireni systemu o integrator
W1=1;
V=(s+10)/s; %presun nuloveho polu a polu -2 na 2x -1
%Korekce pozadavku na citl.funkci filtrem W1
%W1=(s+1)^3/(s*(s+50)^2);
W1=1
%W1=(s+5)/(s+0.2);
sigma(1/(V*W1)); %chovani Sjw pro nizke frekvence
W2=0.01; %mala penalizace rizeni 
sigma(1/(V*W2)); %chovani Ujw pro vysoke frekvence
sigma(G/(V*W2)); %chovani Tjw pro vusoke f
P=ss([W1*V W1*Gext;0 W2;-V -Gext]);
P=minreal(P);
P=mktito(P,1,1);
[K,CL,GAM,INFO] = hinfsyn(P,'display','on');
%kompenzator pro puvodni problem
K0=K/s;
L=G*K0; S=inv(1+L); T=1-S;PS=G*S;CS=K0*S;
close all;
sigma(S,'g',T,'r',GAM/(V*W1),'g-.',GAM*G/ss(V*W2),'r-.',PS,CS)
CL=L/(1+L)
zpk(CL)
pole(CL)
zpk(K0)
%Regulator ma strukturu notch-filter+integrator, tlumeni kolem kmitave
%frekvence by slo zvetsit vahovou funkci velkou kolem tehle frekvence

%Spicka vstupni citlivostni funkce vstupu PS
%Chtelo by to vetsi VF roll-off









%Kmitavy system * zvyseni tlumeni kolem kriticke frekvence 

s=zpk('s');
G=1/(s^2+2*0.1*s+1); 

Gext=G/s; %rozsireni systemu o integrator
W1=1;
%V=(s+10)*(s^2+3*s+1)/(s*(s^2+2*0.1*s+1)); %presun nuloveho polu a polu -2 na 2x -1
V=(s+10)/s;


%Notch filter
N=(s^2+2*0.05*s+1)/(s^2+2*s+1);
W1=1000000/N; %penalizace citlivostni funkce na kritickyzch frekvencich kolem rezonance

%W1=(s+5)/(s+0.2);
%PRI MALE PENALIZACI RIZENI KONSTANTOU W2 BUDOU DVA POLY V -10 DOMINANTNI
sigma(1/(V*W1)); %chovani Sjw pro nizke frekvence
W2=0.00000001; %mala penalizace rizeni 
sigma(1/(V*W2)); %chovani Ujw pro vysoke frekvence
sigma(G/(V*W2)); %chovani Tjw pro vusoke f
P=ss([W1*V W1*Gext;0 W2;-V -Gext]);
P=minreal(P);
P=mktito(P,1,1);
[K,CL,GAM,INFO] = hinfsyn(P,'display','on');
%kompenzator pro puvodni problem
K0=K/s;
L=G*K0; S=inv(1+L); T=1-S;PS=G*S;CS=K0*S;
close all;
sigma(S,'g',T,'r',GAM/(V*W1),'g-.',GAM*G/ss(V*W2),'r-.',PS,CS)
CL=L/(1+L)
zpk(CL)
pole(CL)
zpk(K0)

%Chybi VF roll-off




%Kmitavy system - spravny navrh s aktivnim tlumenim
s=zpk('s');
G=1/(s^2+2*0.1*s+1); 
Gext=G/s; %rozsireni systemu o integrator
W1=1;
V=(s+5)*(s^2+2*s+1)/(s*(s^2+2*0.1*s+1)); %presun zmena kmitavych polu ot.smycky na tlumene v uzavrene
%V=(s+10)/s;
%Notch filter
N=(s^2+2*0.1*s+1)/(s^2+3*s+1);
W1=1; 
%W1=(s+5)/(s+0.2);
sigma(1/(V*W1)); %chovani Sjw pro nizke frekvence
W2=0.001*1000/15*(s+15)/(s+1000); %mala penalizace rizeni 
sigma(W2);
sigma(1/(V*W2)); %chovani Ujw pro vysoke frekvence
sigma(G/(V*W2)); %chovani Tjw pro vusoke f
P=ss([W1*V W1*Gext;0 W2;-V -Gext]);
P=minreal(P);
P=mktito(P,1,1);
[K,CL,GAM,INFO] = hinfsyn(P,'display','on');
%kompenzator pro puvodni problem
K0=K/s;
L=G*K0;
S=inv(1+L); T=1-S;PS=G*S;CS=K0*S;
close all;
sigma(S,'g',T,'r',GAM/(V*W1),'g-.',GAM*G/ss(V*W2),'r-.',PS,CS,'m',GAM/(V*W2),'m-.')
CL=L/(1+L)
zpk(CL)
pole(CL)
zpk(K0)




%Kmitavy system se slabe tlumenym modem, pozadavek na NIZSI sirku pasma
s=zpk('s');
G=1/(s^2+2*0.1*s+1); 
Gext=G/s; %rozsireni systemu o integrator
W1=1;
V=(s+1)^2*(s+0.1)/(s*(s^2+2*0.1*s+1)); %zatlumeni kmitavych polu a + dominantni pol v 0.1
%Korekce pozadavku na citl.funkci filtrem W1
%W1=(s+1)^3/(s*(s+50)^2);
W1=1
%W1=(s+5)/(s+0.2);
sigma(1/(V*W1)); %chovani Sjw pro nizke frekvence
W2=100000*(s+1)/(s+1000); %mala penalizace rizeni 
sigma(1/(V*W2)); %chovani Ujw pro vysoke frekvence
sigma(G/(V*W2)); %chovani Tjw pro vusoke f
P=ss([W1*V W1*Gext;0 W2;-V -Gext]);
P=minreal(P);
P=mktito(P,1,1);
[K,CL,GAM,INFO] = hinfsyn(P,'display','on');
%kompenzator pro puvodni problem
K0=K/s;
L=G*K0; S=inv(1+L); T=1-S;PS=G*S;CS=K0*S;
close all;
sigma(S,'g',T,'r',GAM/(V*W1),'g-.',GAM*G/ss(V*W2),'r-.',PS,'c',CS,'m')
CL=L/(1+L)
zpk(CL)
pole(CL)
zpk(K0)

%Vibrace jsou tlumeny pouze pasivne, dle ocekavani




%DVOJITY INTEGRATOR
s=zpk('s');
G=1/(s^2); 
W1=1;
%V=(s^2+2*sqrt(2)/2*s+1)/s^2; %presun nulovych polu
V=((s+1)^2)/s^2; %presun nulovych polu
%Korekce pozadavku na citl.funkci filtrem W1
%W1=(s+1)^3/(s*(s+50)^2);
W1=1
%W1=(s+5)/(s+0.2);
%PRI MALE PENALIZACI RIZENI KONSTANTOU W2 BUDOU DVA POLY V -10 DOMINANTNI
sigma(1/(V*W1)); %chovani Sjw pro nizke frekvence
W2=(10000/10)*(s+10)/(s+10000); %VF roll-off
sigma(1/(V*W2)); %chovani Ujw pro vysoke frekvence
sigma(G/(V*W2)); %chovani Tjw pro vusoke f
P=ss([W1*V W1*G;0 W2;-V -G]);
P=minreal(P);
P=mktito(P,1,1);
[K,CL,GAM,INFO] = hinfsyn(P,'display','on');
%kompenzator pro puvodni problem
K0=K;
L=G*K0; S=inv(1+L); T=1-S;
close all;
sigma(S,'g',T,'r',GAM/(V*W1),'g-.',GAM*G/ss(V*W2),'r-.')
CL=L/(1+L)
zpk(CL)
pole(CL)
zpk(K0)
%Vychazi derivacni clanek+LP filtr






%Tvorba matice otevrene smycky pro hinfsyn
s=zpk('s');
G=(5)/((s-1)); 
MS=1.05; %max.peak citlivostni funkce
AS=0.1; %utlum poruch na nizkych frekvencich
WS=40; %sirka pasma
W1=(s/MS+WS)/(s+AS*WS);
MT=1.05; %max peak kompl. citl. fce
AT=.005; %utlum na vysokych frekvencich
WT=80; %sirka pasma
W3=(s+WT/MT)/(AT*s+WT);
P=ss([W1 -W1*G;0 W3*G;1 -G]);
P=mktito(P,1,1) %vyssi rad regulatoru, P neni v min. realizaci
P2 = augw(G,W1,[],W3)
P3=minreal(P); %VYHODNE UDELAT MINIMALNI REALIZACI SYSTEMU, JINAK SE ZBYTECNE ZVYSUJE RAD REGULATORU

%[K1,CL,GAM,INFO] = hinfsyn(P);
%[K2,CL,GAM,INFO] = hinfsyn(P2);
[K3,CL,GAM,INFO] = hinfsyn(P3);
L=G*K3; S=inv(1+L); T=1-S;
figure;
sigma(S,'g',T,'r',GAM/W3,'r-.',GAM/W1,'g-.')
figure;
sigma(L,W1,1/W3)






% PRIKLAD KWAKERNAAK - MODIFIKOVANY MIXSYN pro dvojity integrator -
% castecne prirazeni dvojice dominantnich polu
s=zpk('s');
G=1/s^2; 
V=(s^2+sqrt(2)*s+1)/(s^2);
W1=ss(1);

W2=ss(0.1);
P=[W1*V W1*G;0 W2;-V -G];
P=mktito(P,1,1);
P=mktito([W1*V W1*G;0 W2;-V -G],1,1);
P=minreal(P);
[K,CL,GAM,INFO] = hinfsyn(P);
L=G*K; S=inv(1+L); T=1-S;
close all;
sigma(S,'g',T,'r',GAM/W3,'r-.',GAM/W1,'g-.')
CL=L/(1+L);


%PRIKLAD KWAKERNAAK - SYSTEM 1.RADU
%sirka pasma 1 rad/s, regulator s integratorem a VF roll-offem
s=zpk('s');
p=0.1;
G=p/(s+p); 
c=1;
r=0.1;
G0=G/(s*c*(1+r*s))%rozsireny system o integrator
V=(s^2+sqrt(2)*s+1)/(s*(s+p));
W1=1;
W2=1;
P=[W1*V W1*G0;0 W2;-V -G0];
P=ss(P);
P=minreal(P); %BEZ VYPOCTU MINIMALNI REALIZACE DOJDE K PORUSENI PODMINEK NA HODNOSTI MATIC STAVOVE REPRE NUTNE PRO SPRAVNE RESENI
P=mktito(P,1,1);

[K,CL,GAM,INFO] = hinfsyn(P);



%pokus o maticovy zapis rucne nalezene min.realizace ve stavove
%reprezentaci
A=[0 0 p/r;1 -p 0;0 0 -1/r];
B1=[1;sqrt(2)-p;0];
B2=[0;0;1/c];
C1=[0 1 0;0 0 0];
D11=[1;0];
D12=[0;1];
C2=[0 -1 0];
D21=-1;
D22=0;

P2=ss(A,[B1 B2],[C1;C2],[D11 D12;D21 0]);
P2=mktito(P2,1,1)
[K,CL,GAM,INFO] = hinfsyn(P2);

%regulator pro puvodni system
Kp=K/(s*c*(1+r*s));

L=G*Kp; S=inv(1+L); T=1-S;
close all;
sigma(S,'g',T,'r')






%MIMO DESIGN


%2x2 TITO prvni rady, relativne dobre rozvazbene
s=zpk('s');
G=[5/(s+1) 1/(5*s+1);1/(5*s+1) 5/(s+1)];
MS=1.05;AS=.03;WS=5;
W1=(s/MS+WS)/(s+AS*WS);
MT=1.05;AT=.05;WT=20;
W3=(s+WT/MT)/(AT*s+WT);
% Compute the H-infinity mixed-sensitivity optimal sontroller K1
[K1,CL1,GAM1]=mixsyn(G,W1,[],W3);
% Next compute and plot the closed-loop system.
% Compute the loop L1, sensitivity S1, and comp sensitivity T1:
L1=G*K1;
I=eye(size(L1));
%S1=feedback(I,L1); 
 S1=inv(I+L1);
T1=I-S1;
CS=K1*S1;
PS=G*S1;
%casove odezvy
step(T1,1.5);
% frekvencni odezvy
figure;
%Citlivostni a kompl.citl. fce
sigma(S1,'r',GAM1/W1,'r--',T1,'g',GAM1/W3,'g--')
figure;
%Citlivostni funkce sumu a vstupu
sigma(CS,PS);

%CHYBI INTEGRATOR A VF ROLL-OFF

%Zopakovani predchoziho prikladu pres hinfsyn, pro ujasneni dimenzi matic
s=zpk('s');
G=[5/(s+1) 1/(5*s+1);1/(5*s+1) 5/(s+1)];
MS=1.05;AS=.03;WS=5;
W1=(s/MS+WS)/(s+AS*WS);
MT=1.05;AT=.05;WT=20;
W3=(s+WT/MT)/(AT*s+WT);
%Korekce dimenzi pro MIMO pripad 
W1=W1*eye(2);
W3=W3*eye(2);

P=ss([W1 -W1*G;zeros(2,2) W3*G;eye(2) -G]);
P=minreal(P);
P=mktito(P,2,2) %vyssi rad regulatoru, P neni v min. realizaci
%P2 = AUGW(G,W1,[],W3)
%P3=minreal(P); %VYHODNE UDELAT MINIMALNI REALIZACI SYSTEMU, JINAK SE ZBYTECNE ZVYSUJE RAD REGULATORU
[K1,CL,GAM,INFO] = hinfsyn(P);
%[K1,CL1,GAM1]=mixsyn(G,W1,[],W3);
% Next compute and plot the closed-loop system.
% Compute the loop L1, sensitivity S1, and comp sensitivity T1:
L1=G*K1;
I=eye(size(L1));
%S1=feedback(I,L1); 
 S1=inv(I+L1);
T1=I-S1;
CS=K1*S1;
PS=G*S1;

close all;
%casove odezvy
step(T1,1.5);
% frekvencni odezvy
figure;
%Citlivostni a kompl.citl. fce
sigma(S1,'r',GAM1/W1,'r--',T1,'g',GAM1/W3,'g--')
figure;
%Citlivostni funkce sumu a vstupu
sigma(CS,PS);

%Stejny vysledek jako s prikazem mixsyn



%ODLISNE VAZENI JEDNOTLIVYCH KANALU
s=zpk('s');
G=[5/(s+1) 1/(5*s+1);1/(5*s+1) 5/(s+1)];
MS=1.05;AS=.03;
WS1=5; %sirka pasma pro 1.kanal
WS2=50; %2.kanal
W11=(s/MS+WS)/(s+AS*WS);
W12=(s/MS+WS2)/(s+AS*WS2);
MT=1.05;AT=.05;
WT1=20;
WT2=80;
W31=(s+WT1/MT)/(AT*s+WT1);
W32=(s+WT2/MT)/(AT*s+WT2);
%Korekce dimenzi pro MIMO pripad 
W1=[W11 0;0 W12];
W3=[W31 0;0 W32];
P=ss([W1 -W1*G;zeros(2,2) W3*G;eye(2) -G]);
P=minreal(P);
P=mktito(P,2,2) %vyssi rad regulatoru, P neni v min. realizaci
%P2 = AUGW(G,W1,[],W3)
%P3=minreal(P); %VYHODNE UDELAT MINIMALNI REALIZACI SYSTEMU, JINAK SE ZBYTECNE ZVYSUJE RAD REGULATORU
[K1,CL,GAM,INFO] = hinfsyn(P);
%[K1,CL1,GAM1]=mixsyn(G,W1,[],W3);
% Next compute and plot the closed-loop system.
% Compute the loop L1, sensitivity S1, and comp sensitivity T1:
L1=G*K1;
I=eye(size(L1));
%S1=feedback(I,L1); 
 S1=inv(I+L1);
T1=I-S1;
CS=K1*S1;
PS=G*S1;
close all;
%casove odezvy
step(T1,1.5);
% frekvencni odezvy
figure;
%Citlivostni a kompl.citl. fce prvniho kanalu
sigma(S1(1,1),'r',GAM1/W11,'r--',T1(1,1),'g',GAM1/W31,'g--');
figure;
%Druheho
sigma(S1(2,2),'m',GAM1/W12,'m--',T1(2,2),'c',GAM1/W32,'c--')
%Mimodiagonalni citlivostni funkce
figure;
sigma(S1(1,2),S1(2,1),T1(1,2),T1(2,1));

figure;
%Citlivostni funkce sumu a vstupu
sigma(CS,PS);

%Spravne rozvazbeni smycek, poradl ale chybi integrator a VF roll-off




%PRIDANI INTEGRATORU
s=zpk('s');
G=[5/(s+1) 1/(5*s+1);1/(5*s+1) 5/(s+1)];
Gext=G/s; %rozsireni o integrator
%Gext=G.*[1/s 1;1 1/s];
V=(s+10)/s; %presun nuloveho polu do -10 v obou kanalech
W1=1;
sigma(1/(V*W1)); %omezeni citl.fce pro oba kanaly, sirka pasma 10 rad/s
W2=0.02;
sigma(1/(V*W2));
%Korekce dimenzi pro 2x2 MIMO, penalizace diagonalnich prvku
W1=W1*eye(2);
W2=W2*[1 0;0 1];
P=ss([W1*V W1*Gext;zeros(2,2) W2;-V*eye(2) -Gext]);
P=minreal(P);
P=mktito(P,2,2);
[K,CL,GAM,INFO] = hinfsyn(P,'display','on');
%Kompenzator pro puvodni system bez integratoru
K0=K/s;
L1=G*K0;
I=eye(size(L1));
%S1=feedback(I,L1); 
 S1=inv(I+L1);
T1=I-S1;
CS=K0*S1;
PS=G*S1;
close all;
%casove odezvy
step(T1,1.5);
% frekvencni odezvy
figure;
%Citlivostni a kompl.citl. fce
sigma(S1,'r',GAM/(V*W1(1,1)),'r--',T1,'g',GAM/(V*W2(1,1)),'g--');
%Mimodiagonalni citlivostni funkce
figure;
sigma(S1(1,2),S1(2,1),T1(1,2),T1(2,1));
figure;
%Citlivostni funkce sumu a vstupu
sigma(CS,PS);









%VF rolloff a integrator
s=zpk('s');
G=[5/(s+1) 1/(5*s+1);1/(5*s+1) 5/(s+1)];
Gext=G/s; %rozsireni o integrator
%Gext=G.*[1/s 1;1 1/s];
V=(s+10)/s; %presun nuloveho polu do -10 v obou kanalech
W1=1;
sigma(1/(V*W1)); %omezeni citl.fce pro oba kanaly, sirka pasma 10 rad/s
%W2=0.000001*(s+15)/(s+1000); %penalizace rizeni pro omezeni sirky pasma a VF roll-off
W2=0.05*100/15*(s+15)/(s+100);
%W2=0.1;
sigma(1/(V*W2));
%Korekce dimenzi pro 2x2 MIMO, penalizace diagonalnich prvku
W1=W1*eye(2);
W2=W2*[1 0;0 1];
P=ss([W1*V W1*Gext;zeros(2,2) W2;-V*eye(2) -Gext]);
P=minreal(P);
P=mktito(P,2,2);
[K,CL,GAM,INFO] = hinfsyn(P,'display','on');
%Kompenzator pro puvodni system bez integratoru
K0=K/s;
L1=G*K0;
I=eye(size(L1));
%S1=feedback(I,L1); 
 S1=inv(I+L1);
T1=I-S1;
CS=K0*S1;
PS=G*S1;
close all;
%casove odezvy
step(T1,1.5);
% frekvencni odezvy
figure;
%Citlivostni a kompl.citl. fce
sigma(S1,'r',GAM/(V*W1(1,1)),'r--',T1,'g',GAM/(V*W2(1,1)),'g--');
%Mimodiagonalni citlivostni funkce
figure;
sigma(S1(1,2),S1(2,1),T1(1,2),T1(2,1));
figure;
%Citlivostni funkce sumu a vstupu
sigma(CS,PS);


%REGULATOR MA INTEGRATOR A DODATECNY VF ROLL OFF






























%NESTABILNI SYSTEM - 1 nestabilni pol

%SPATNY NAVRH - REGULATOR KRATI NESTABILNI POL NESTABILNI NULOU
%Pri minimalni chybe v modelu systemu bude uzavr. smycka nestabilni !!!

s=zpk('s');
G=1/((s-1)*(s+1)^3); 
Gext=G/s; %rozsireni systemu o integrator
W1=1;
%V=(s^2+2*s+1)/(s*(s+2)); %presun nuloveho polu a polu -2 na 2x -1
V=(s^2+2*1*s+1)/(s*(s-1));
%Korekce pozadavku na citl.funkci filtrem W1
%W1=(s+1)^3/(s*(s+50)^2);
W1=1;
%PRI MALE PENALIZACI RIZENI KONSTANTOU W2 BUDOU DVA POLY V -10 DOMINANTNI
sigma(1/(V*W1)); %chovani Sjw pro nizke frekvence
W2=0.01; %mala penalizace rizeni 
P=ss([W1*V W1*Gext;0 W2;-V -Gext]);
P=minreal(P);
P=mktito(P,1,1);
[K,CL,GAM,INFO] = hinfsyn(P,'display','on');
%kompenzator pro puvodni problem
K0=K/s;
L=G*K0; S=inv(1+L); T=1-S;
close all;
sigma(S,'g',T,'r',GAM/(V*W1),'g-.',GAM*G/ss(W2),'r-.')
CL=L/(1+L)
zpk(CL)
pole(CL)
zpk(K0)