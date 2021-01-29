function [han,modele]=ident_wynik(czas,mrac_ucz,wyjscia,sterowanie,referencje,varargin)
legU={'$k_r$','$k_x$'};
legT={'MRAC','PID'};
xL='Numer cyklu';
yL='Czêstotliwoœæ log_{10}[rad/s]';
zP='[\circ]';
zM='[dB]';
han(1)=figure;
N=2;

plot(czas,mrac_ucz(:,1))
hold on
plot(czas,mrac_ucz(:,2))
if size(mrac_ucz,2)<3
estetyka(legU,'Proces uczenia')
else
plot(czas,mrac_ucz(:,3))
plot(czas,mrac_ucz(:,4))
estetyka(compose('$\\theta_%d$',1:4),'Proces uczenia')
end
% 
% subplot(1,2,2)
% plot(czas,mracN_ucz(:,1))
% hold on
% plot(czas,mracN_ucz(:,2))
% legend(legU,'location','best')
% title('mracN')
% hold off


han(2)=figure;
subplot(2,1,1)
hold on
for id=1:N
plot(czas,wyjscia(:,id))
end

snr=mag2db(rms(referencje(:,1))/rms(referencje(:,1)-referencje(:,2)));

estetyka([],sprintf('y(t), SNR pomiarowy: %2.2f [dB]',snr))


subplot(2,1,2)
hold on
for id=1:N
plot(czas,sterowanie(:,id))
end
estetyka(legT,'u(t)')
%}
% uchyby=referencje(:,2)-wyjscia; %2-bez szumu
% 
% subplot(3,1,3)
% hold on
% for id=1:N
% plot(czas,uchyby(:,id))
% end
% estetyka(legT,'e(t)')

%bode estymacja, odpowiedz skokowa
if nargin>5
    par=varargin{1};
    
    Ts=par.Ts;
    rzadl=par.rzadl;
    rzadm=par.rzadm;
    opozn=par.opozn;
    znacznik=par.znacznik;
    len=length(znacznik)-1;

    modele=cell(len,2);
    omeg=logspace(-1,log10(.5/Ts*2*pi),100);
    modele{1,1}.omega=omeg;
    
    %lepszy parfor
    for id=1:len
       bity = czas>znacznik(id) & czas<znacznik(id+1);
       segCz =  czas(bity);
       t = segCz(1):Ts:segCz(end);
       
       segX = referencje(bity,1);
       segXi=interp1(segCz,segX,t,'spline');
       %plot(mag2db(abs(fft(segYi)./fft(segXi)))) %szumy, nieczytelne

       for idj=1:N
       segY = wyjscia(bity,idj);
       segYi=interp1(segCz,segY,t,'spline');
       
       dane=iddata(segYi',segXi',Ts);
       dane.InterSample='foh';
       
              
       model=tfest(dane,rzadm,rzadl,opozn);
       [~,fit]=compare(model,dane);
       if fit<50
%             if idj==1
%            fprintf('\nmrac c: %d, l: %d  m: %d\n',id,rzadl,rzadm)
%             else
%            fprintf('\npid c: %d, l: %d  m: %d\n',id,rzadl,rzadm)
%             end
       end

       zajety=fit<80;
       
       while zajety

    %{
       
       if fit<90 || idj==1
       figure
       if idj==1
           title(sprintf('\nmrac c: %d, l: %d  m: %d\n',id,rzadl,rzadm))
       else
           title(sprintf('\npid c: %d, l: %d  m: %d\n',id,rzadl,rzadm))
       end
%        compare(model,dane);
%        pause
%            zajety=1;
%            opozn=NaN;
%        else

%            zajety=0;
       end
      %}

%           figure(100)
%            compare(model,dane)

%        if idj==1
%            title(sprintf('mrac c: %d, l: %d  m: %d',id,rzadl,rzadm))
%        else
%            title(sprintf('pid c: %d, l: %d  m: %d',id,rzadl,rzadm))
%        end
%        zajety=input('ok? [0]\npoprawa? [1]\n>>');
%            rzadm=input('bieguny: ');
%            rzadl=input('zera: ');
           if fit<55 %dla pid po rozregulowaniu
           model=tfest(dane,rzadm+1,rzadl+1,opozn); %3,1
%            pause
%            compare(model,dane)
           end
%            sprintf('%d c: %d, l: %d  m: %d',idj,id,rzadl,rzadm)
%            pause
           zajety=0;


%        rzadl=par.rzadl;
%        rzadm=par.rzadm;
       end
       [mag,ph,~]=bode(model,omeg);
       
       modele{id,idj}.mag=mag2db(mag);
%        modele{id,idj}.ph=ph;

       modele{id,idj}.ph=mod(ph+180,360)-180;
       
       modele{id,idj}.estym=model;
       end
    end%for id
    
    %{
    for idj=1:N
    [DX,DY]=meshgrid(1:len,log10(omeg));
    
    for id=len:-1:1
        magS(:,id)=modele{id,idj}.mag;
        phaS(:,id)=modele{id,idj}.ph;
    end
    
    han(4+2*idj-1)=figure;
    surf(DX,DY,magS)
    view([120 45])
    title(['Charakterystyki amplitudowe przebiegów ' legT{idj}])
    xlabel(xL);    ylabel(yL);    zlabel(zM);
    han(4+2*idj)=figure;
    surf(DX,DY,phaS)
    view([120 45])
    title(['Charakterystyki fazowe przebiegów ' legT{idj}])
    xlabel(xL);    ylabel(yL);    zlabel(zP);
    end%for idj
    %}
elseif nargout>1
    modele=[];   
end%nargin




function estetyka(leg,tyt)
    if ~isempty(leg)
    legend(leg,'location','best','Interpreter','latex')
    end
    title(tyt)
    grid on
    grid minor
    xlabel('t [s]','Interpreter','latex')
    hold off
end
end