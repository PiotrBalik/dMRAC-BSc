function [wskazniki]=odp_czas(czas,wyjscia,referencje,varargin)

par=varargin{1};

Ts=par.Ts;
znacznik=par.znacznik;
N=size(wyjscia,2);

%wskazniki regulacji
len=length(znacznik)-1;
wsk=zeros(4,len,2);


for id=1:len
    bity = czas>znacznik(id) & czas<znacznik(id+1);
    segCz =  czas(bity);
    t = segCz(1):Ts:segCz(end);
    len=length(t);
    
    ref = referencje(bity,2);
    zmiana=0;
    if ref(round(len/2))<0
        ref=-ref;
        zmiana=1;
    end
    ref=interp1(segCz,ref,t,'spline');
    
    
    for idj=1:N
        segY = wyjscia(bity,idj);
        if zmiana
            segY=-segY;
        end
        y=interp1(segCz,segY,t,'spline');
        
        %suma e^2
        e=y-ref;
        wsk(1,id,idj)=sum(e.^2)/len;
        
        
        DC=mean(ref(round(len*0.95):len)); %ostatnie 5% probek
        
        %narastanie
        wzr1=find(y>DC*0.1,1,'First');
        wzr2=find(y>DC*0.9,1,'First');
        wsk(2,id,idj)=t(wzr2)-t(wzr1);
        
        %ustalanie dla 9%
        ust=find(abs(e)>DC*0.09,1,'Last');
        wsk(3,id,idj)=t(ust)-t(1);
        
        %przereg
        wsk(4,id,idj)=(max(y)-DC)/DC*100;
%         if id>7
%         plot(t,y,t,ref,t,abs(e),t,ones(len,1)*DC*0.08)
%         title([sprintf('ster %d,cykl %d ',idj,id) sprintf('p=%2.2f, ',wsk(:,id,idj))])
%         pause
%         end
        
    end%for idj
end%for id

wskazniki=wsk;
end