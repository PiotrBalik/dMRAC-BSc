%% inicjalizacja
% set_param('mainLinR1/ref/resztaVramp','Value','1') 
% set_param('mainLinR1/ref/sinVpuls','Value','1') %sinus

% set_param('mainLin/ref/sinVpuls','Value','-1') %puls
% set_param('mainLin','StopTime','200')

% set_param('mainLin','IgnoredZcDiagnostic','warning')
% set_param('mainLin','MaskedZcDiagnostic','warning');

% set_param('mainLinR1','IgnoredZcDiagnostic','none');
% set_param('mainLinR1','MaskedZcDiagnostic','none');

% set_param('mainNLin','IgnoredZcDiagnostic','none');
% set_param('mainNLin','MaskedZcDiagnostic','none');

clear variables
close all

% %parpool('local',2)
% a=parcluster('local');
%parpool('local',a.NumWorkers)
Tsim=120;

han=[];
legB={'pid','mrac','obiekt uk³. otwarty'};
open_system('mainLinR1.slx')
open_system('mainLinR2.slx')

load('nastawy.mat')
%analiza segmentow
parametry.Ts=0.01;
okres=10;
parametry.znacznik=0:okres:Tsim;
parametry.opozn=0;
% Riccati
%A'X + XA - XBB'X + Q = 0
% do P*A+A'P + Q = 0
% B=zeros(n,1);
% P=care(NAref,B,Q);

wybor=input('[1] - liniowy\n[2] - liniowy\n[3] - nieliniowy\n>>','s');
%% liniowy model 1-rzad

if any(wybor=='1')
sys1=zpk([],-1/.34,.9);
sys1=ss(sys1);

Am=sys1.A;  Bm=sys1.B;  Cm=sys1.C;
xmodel0=0;

%do estymacji
parametry.rzadm=1;
parametry.rzadl=0;

N=1/C1.Tf;
D=C1.Kd;
I=C1.Ki;
PPid=C1.Kp;
% PPid=5.56;
% I=16;
% D=1.16;
% N=8e-04;
s=tf('s');
pid=PPid+I/s+D*N/(1+N/s);


%martwa strefa
epsilon=0.1;
szum=3e-3;



% 2-gi rzad, relatywny 1
Aref=-5;
mrac=zpk([],Aref,1);
% mrac=zpk(Aref,[Aref+[1i -1i]*3.5],1); %dodatkowe zero
mrac=ss(mrac/dcgain(mrac));

NAref=mrac.A;   NBref=mrac.B;   NCref=mrac.C;
Xref0=0;
Kx0=0;
Kr0=0;
p0=1;
lambdaR=10;
lambdaX=10;

pid_z=feedback(pid*sys1,1);

% bode(pid_z,mrac,sys1);
% legend(legB)
% wykres=gcf;
% han=[han wykres];


Auszk=-Am+2; %uszkodzenie
Buszk=Bm*.55;
TuszkA=okres*4;     TuszkB=okres*8;
Azaburz=[];     Bzaburz=[];
Azaburz.time=[1:Tsim]';
Azaburz.signals.values=[(Azaburz.time>=TuszkA)*Auszk]; %skokowe

Bzaburz.time=[1:Tsim]';
Bzaburz.signals.values=[-(Bzaburz.time>=(TuszkB))*Buszk]; %skokowe


set_param('mainLinR1/ref/resztaVramp','Value','1') 
set_param('mainLinR1/ref/sinVpuls','Value','-1') %puls
sim('mainLinR1.slx')

[han1, modele1s]=ident_wynik(tout,mrac_ucz,wyjscia,sterowanie,referencje,parametry);
wsk1s=odp_czas(tout,wyjscia,referencje,parametry);
han=[han han1];

%wizualizacja modeli
odniesienia1={mrac,pid_z,sys1};
cykle=[1 TuszkA/okres+1 TuszkB/okres+1];
han=[han wykresyBode(modele1s,odniesienia1,cykle)];


%narastajace
Azaburz.signals.values=[zeros(1,TuszkA) linspace(0,Auszk,TuszkB-TuszkA) Auszk*ones(1,Tsim-TuszkB)]'; %narastajace
Bzaburz.signals.values=[zeros(1,TuszkB) linspace(0,-Buszk,Tsim-TuszkB)]';


sim('mainLinR1.slx')

[han1, modele1r]=ident_wynik(tout,mrac_ucz,wyjscia,sterowanie,referencje,parametry);
wsk1r=odp_czas(tout,wyjscia,referencje,parametry);
han=[han han1];

cykle=[1 TuszkB/okres+1 Tsim/okres];
han=[han wykresyBode(modele1r,odniesienia1,cykle)];

end

%% liniowy model 2-rzad
if any(wybor=='2')
% silnik DC
% J = 0.01;
% b = 0.1;
% K = 0.01;
% R = 1;
% L = 0.5;
% Am = [-b/J   K/J
%     -K/L   -R/L];
% Bm = [0
%     1/L];
% Cm = [1   0];

% zawieszenie
m=5; %[kg]
c=100; %[N*s/m]
k=40; %[N/m]
Am=[-c/m -k/m;
    1 0];
Bm=[1/m;
    0];
Cm=[0 1];
sys2=ss(Am,Bm,Cm,[]);

xmodel0=[0 0]';
szum=2e-4;      epsilon2=0.02;

%do estymacji
parametry.rzadm=2;      parametry.rzadl=0;

N=1/C2.Tf;
D=C2.Kd;
I=C2.Ki;
PPid=C2.Kp;
% PPid=106.2881;
% I=53.3851;
% D=22.7393;
% N=2.0470e-05;
pid=PPid+I/s+D*N/(1+N/s);
pid_z=feedback(pid*sys2,1);

Aref=-0.4;
mrac=zpk([],[Aref-2+[1i -1i]*Aref/1.5],1);
mrac=ss(mrac/dcgain(mrac));

NAref=mrac.A;   NBref=mrac.B;   NCref=mrac.C;

Nxref0=zeros(2,1);
w10=0;  w20=0;
phi0=zeros(4,1);    p0=1;
Th0=zeros(4,1); %Th0(3)=-60;
F=-2;   g=.5;
Gamma=10*eye(4);


% bode(pid_z,mrac,sys2);
% legend(legB)
% wykres=gcf;
% han=[han wykres];

Auszk=[0 k/m; 0 0];
Buszk=-Bm/2;
TuszkA=okres*4;     TuszkB=okres*8;
tvec=[1:Tsim];

Azaburz=[];     Bzaburz=[];
Azaburz.time=tvec';

%skokowe
tmp=repmat(Auszk,[1 1 length(tvec)]);
tmp(:,:,(tvec<TuszkA))=repmat(zeros(2),[1 1 sum(tvec<TuszkA)]);

Azaburz.signals.values=tmp;
Azaburz.signals.dimensions=[2 2];

Bzaburz.time=tvec';
tmp=repmat(Buszk,[1 1 length(tvec)]);
tmp(:,:,(tvec<TuszkB))=repmat(zeros(2,1),[1 1 sum(tvec<TuszkB)]);

Bzaburz.signals.values=tmp;
Bzaburz.signals.dimensions=[2 1];


set_param('mainLinR2/ref/resztaVramp','Value','1') 
set_param('mainLinR2/ref/sinVpuls','Value','-1') %puls
sim('mainLinR2.slx')

[han5, modele2s]=ident_wynik(tout,mrac2_ucz,wyjscia,sterowanie,referencje,parametry);
wsk2s=odp_czas(tout,wyjscia,referencje,parametry);
han=[han han5];

%wizualizacja modeli
odniesienia2={mrac,pid_z,sys2};
cykle=[TuszkA/okres TuszkA/okres+1 TuszkB/okres+1];
han=[han wykresyBode(modele2s,odniesienia2,cykle)];

NN=TuszkB-TuszkA;
%narastajace
tmp=repmat(Auszk,[1 1 Tsim]);
tmp(:,:,(tvec<TuszkA))=repmat(zeros(2),[1 1 TuszkA-1]);
tmp2=reshape(tmp(:,:,(TuszkA<=tvec & tvec<TuszkB)),[4 NN]);
tmp2=tmp2.*linspace(0,1,NN);
tmp(:,:,(TuszkA<=tvec & tvec<TuszkB))=reshape(tmp2,[2 2 NN]);

Azaburz.signals.values=tmp;

NN=Tsim-TuszkB+1;
tmp=repmat(Buszk,[1 1 Tsim]);
tmp(:,:,(tvec<TuszkB))=repmat(zeros(2,1),[1 1 TuszkB-1]);
tmp2=reshape(tmp(:,:,(tvec>=TuszkB)),[2 NN]);
tmp2=tmp2.*linspace(0,1,NN);
tmp(:,:,(TuszkB<=tvec))=reshape(tmp2,[2 1 NN]);

Bzaburz.signals.values=tmp;

sim('mainLinR2.slx')

[han6, modele2r]=ident_wynik(tout,mrac2_ucz,wyjscia,sterowanie,referencje,parametry);
wsk2r=odp_czas(tout,wyjscia,referencje,parametry);
han=[han han6];

%wizualizacja modeli
cykle=[1 TuszkB/okres+1 Tsim/okres];
han=[han wykresyBode(modele2r,odniesienia2,cykle)];

end
%% nieliniowy
if any(wybor=='3')
nln=0;
if nln==1
%M. Boufadene, Nonlinear Control Systems Using MATLAB, CRC Press 2019
%dx = [x2; -sin(x1) + u] --wahadlo
% pozycja, predkosc
xmodel0=[pi/6 -0.1]';
fmodel=@(x) [x(2) -sin(x(1))+x(3)]'; % ostatnie stany to wejscia
% wmodel=@(x) [(mod(x(1)+pi,2*pi)-pi)*180/pi x(2)]';
wmodel=@(x) (mod(x(1)+pi,2*pi)-pi)*180/pi;
elseif nln==2
%Vaidyanathan, Sundarapandian & Member, Non. (2007). Output Regulation of Van der Pol Oscillator. 88. 
%dx = [x2; -x1 -mu x2(x1^2-1) + u]
mu=-1;
xmodel0=[2 0.1]';
fmodel = @(x) [x(2); -x(1)-mu*x(2)*(x(1)^2-1)+x(3)];
% fmodel = @(x) [x(2); -x(1)-mu*x(2)+x(3)];

wmodel=@(x) [x(1)];
% wmodel=@(x) [x(1) x(2)]';
end
m=5; %[kg]
c=100; %[N*s/m]
k=40; %[N/m]
%   predk  przes
Am=[-c/m -k/m;
    1 0];
Bm=[1/m;
    0];
Cm=[0 1];
sys2=ss(Am,Bm,Cm,[]);

AmN=Am;
AmN(1)=0;

%nieloniowowsc tlumika
fmodel=@(x) AmN*x(1:2)+Bm*x(3) + [-c/m*atan(x(1)^3) ;0];
wmodel=@(x) x(2);

xmodel0=[0 0]';
szum=2e-4;      epsilon2=0.02;
Tsim=80;

Auszk=[0 k/m; 0 0];
Buszk=Bm*0;

TuszkA=okres*4;     TuszkB=okres*4;
tvec=[1:Tsim];

Azaburz=[];     Bzaburz=[];
Azaburz.time=tvec';

%skokowe
tmp=repmat(Auszk,[1 1 length(tvec)]);
tmp(:,:,(tvec<TuszkA))=repmat(zeros(2),[1 1 sum(tvec<TuszkA)]);

Azaburz.signals.values=tmp;
Azaburz.signals.dimensions=[2 2];

Bzaburz.time=tvec';
tmp=repmat(Buszk,[1 1 length(tvec)]);
tmp(:,:,(tvec<TuszkB))=repmat(zeros(2,1),[1 1 sum(tvec<TuszkB)]);

Bzaburz.signals.values=tmp;
Bzaburz.signals.dimensions=[2 1];

%do estymacji
parametry.rzadm=2;      parametry.rzadl=0;

N=1/C2.Tf;
D=C2.Kd;
I=C2.Ki;
PPid=C2.Kp;
% PPid=106.2881;
% I=53.3851;
% D=22.7393;
% N=2.0470e-05;
pid=PPid+I/s+D*N/(1+N/s);
pid_z=feedback(pid*sys2,1);

Aref=-0.4;
mrac=zpk([],[Aref-2+[1i -1i]*Aref/1.5],1);
mrac=ss(mrac/dcgain(mrac));

NAref=mrac.A;   NBref=mrac.B;   NCref=mrac.C;

Nxref0=zeros(2,1);
w10=0;  w20=0;
phi0=zeros(4,1);    p0=1;
Th0=zeros(4,1); %Th0(3)=-60;
F=-2;   g=.5;
Gamma=10*eye(4);


% bode(pid_z,mrac,sys2);
% legend(legB)
% wykres=gcf;
% han=[han wykres];

sim('mainNLinR2.slx')


parametry.Ts=0.01;
okres=10;
parametry.znacznik=0:okres:Tsim;
parametry.opozn=0;
[han5, modele2n]=ident_wynik(tout,mrac2_ucz,wyjscia,sterowanie,referencje);
wsk2n=odp_czas(tout,wyjscia,referencje,parametry);
han=[han han5];

latexport(han,{'uczN','odpN'})

    for id=1:4
        %mrac - 1, pid - 2
        tabela(2*(id-1)+1,:)=wsk2n(id,:,1);
        tabela(2*(id-1)+2,:)=wsk2n(id,:,2);

    end
    arr2lat(tabela,'modelN%d');



% epsTh=0.5;
% thMax=2;
% FP=@(th) (norm(th)^2-thMax^2)/(epsTh*thMax^2);
% gradFP=@(th) 2*th/(epsTh*thMax^2);
% projection=@(th,y) y-(gradFP(th)*gradFP(th)')/norm(gradFP(th))^2*y*FP(th)*((FP(th)>0) && (y'*gradFP(th) > 0));
end
%% export
%delete(gcp)
% savefig(han,'PeaksFile.fig','compact')
% load('wyniki_19_01_lepsze.mat')
if all(wybor~='3')
eksport=input('eksport? [1]: ');

if export==1
    %tytu={'bode','k1','y1','u1','e1','sm1','sm2','sp1','sp2','char1','k2','y2','u2','e2','sm1','sm2','sp1','sp2','char2'};
    glown={'ucz','odp','amp','faz'};
    cykle={'1s','1r','2s','2r'};
    tytu=cell(1,16);
    for idi=1:4
        for idj=1:4
            tytu{(idi-1)*4+idj}=[glown{idj} cykle{idi}];
        end
    end
latexport(han,tytu)
styl={'mrac','pid'};
wsks=cell(4,1);
wsks{1}=wsk1s;
wsks{2}=wsk1r;
wsks{3}=wsk2s;
wsks{4}=wsk2r;

tyt_wsk={'inte','tr','ts','ov'};
for it=[1 3]

    for id=1:4
        %mrac - 1, pid - 2
        tabela(4*(id-1)+1,:)=wsks{it}(id,:,1);
        tabela(4*(id-1)+2,:)=wsks{it}(id,:,2);
        
        tabela(4*(id-1)+3,:)=wsks{it+1}(id,:,1);
        tabela(4*(id-1)+4,:)=wsks{it+1}(id,:,2);
    end
%     mat2lat(tabela,sprintf('model%d',ceil(it/2)));
    arr2lat(tabela,sprintf('model%d',ceil(it/2)));

end

% nazwy={   'Twskinte'
%     'Twsktr'
%     'Twskts'
%     'Twskov'
% }
% for id=1:4
    %generowanie tabel
%     polec{id}=sprintf('%s2=array2table(%s,''RowNames'',tyt_tab);',nazwy{id},nazwy{id});
    
    %podpisywanie
%       polec{id}=sprintf('%s2.Properties.VariableNames=NN;',nazwy{id});

%    eval(['wsk1' cykle{id} 'T=table();'   ]);
% end

% fprintf('%s\n',polec{:})

% tyt_tab=cell(1,8);
% tyt_tab(1:2:8)=compose('MRAC%d',1:4);
% tyt_tab(2:2:8)=compose('PID%d',1:4);

%load('Tsk.mat')
% splitvars(table(Twskinte),compose('Cykl: %d',1:12))

%Twsk.Properties.VariableNames = compose('Cykl: %d',1:12);

end
end