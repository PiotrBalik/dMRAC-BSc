function [han]=wykresyBode(modele,odniesienie,cykle)
%odniesienie: 1x3, mrac,pid,system
N=length(modele);
legN={'MRAC','PID','Obiekt'};

omega=modele{1,1}.omega;
omega=log10(omega);
% [ms,phs]=bode(odniesienie{3},10.^omega);

p='Position';
pn=[400 200 900 400];

han(1)=figure(p,pn);
podwykres({'Char. amplitudowe ', 'A(f)\enskip[dB]'},1)

han(2)=figure(p,pn);
podwykres({'Char. fazowe ','\varphi(f)\enskip[deg]'},0)


function podwykres(tekst,typ)
for idj=1:2
[mr,phr]=bode(odniesienie{idj},10.^omega);
[mo,pho]=bode(odniesienie{3},10.^omega);
subplot(1,2,idj)
hold on
% plot(omega,mag2db(ms(:)))
if typ
semilogx(omega,mag2db(mo(:)),'-.')
mr=mr(:);
semilogx(omega(1:2:end),mag2db(mr(1:2:end)),'+','MarkerSize',3)
%adaptacja dla osi...
ylim([round(mag2db(min(mo(:)))*0.75,-1) 15])
else
semilogx(omega,pho(:),'-.')
phr=phr(:);
semilogx(omega(1:2:end),phr(1:2:end),'+','MarkerSize',3)
end

if idj<3
for id=cykle %1:2:N
    subplot(1,2,idj)
    if typ
    semilogx(omega,modele{id,idj}.mag(:))
    else
    semilogx(omega,modele{id,idj}.ph(:))
    end
end
end%idj<3
% legE={'Referencja',legT{cykle}};
% legE={'Uk. otwarty','Referencja',legT{cykle}};
subplot(1,2,idj)
estetyka([ tekst{1} legN{idj} ],tekst{2})


%dla amplitud
if typ
   %legenda
   if idj==1
   tmp=flip(compose('Cykl: %d',cykle));
   tmp{end+1}='Wymagania';
   tmp{end+1}=legN{3}; %obiekt
   legend(flip(tmp),'location','southwest')
   end
   miniwykres(idj)
end


end%for idj
end%funkcja

function miniwykres(idj)
    a=gca;
%     yl=a.YLim(1);
%     ym=a.YLim(2);
    a=a.Children;
    n=length(a);
    
%     offset=round((round(max(a(end-1).YData))-yl)/(ym*1.08-yl),2);
%     annotation('arrow',[0.34 0.22]+(idj-1)*.44,[0.84 offset])

%kradzione, normalizacja wspolrzednych
%https://stackoverflow.com/users/471614/bdecaf
    ax = axis;
    ap = get(gca,'Position');
    xo = [1.6 0.2];
    yo = [ax(4)/2 ax(4)/6];
    xp = (xo-ax(1))/(ax(2)-ax(1))*ap(3)+ap(1);
    yp = (yo-ax(3))/(ax(4)-ax(3))*ap(4)+ap(2);
    annotation('arrow',xp,yp)

    %dwie pozycje
    b=axes('Position',[0.37+(idj-1)*.44 0.725 0.1 0.2],'Layer','top');
    plot(b,rand(2,n))
    
    for linia=1:n %kopiowanie
        b.Children(linia).XData=a(linia).XData;
        b.Children(linia).YData=a(linia).YData;
        b.Children(linia).Color=a(linia).Color;
        b.Children(linia).LineStyle=a(linia).LineStyle;
        b.Children(linia).Marker=a(linia).Marker;
        b.Children(linia).MarkerSize=a(linia).MarkerSize;
    end
    ylim(b,[-3 3])
    xlim(b,[-1 0.7])
    grid(b,'on')
    grid(b,'minor')
end
function estetyka(tyt,yL)
    title(tyt)
    grid('on')
    grid('minor')
    xlabel('$$log_{10}(f)\enskip [rad/s]$$','Interpreter','latex')
    ylabel(['$$' yL '$$'],'Interpreter','latex')
    hold('off')
end
end