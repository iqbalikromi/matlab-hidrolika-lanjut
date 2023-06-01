clear
clf

So=0;
g=9.81;

m=0;
% n(i)=0;
Qo=0;
dt=0.001;
time=60;
tm=time/dt;
tplot=1;
dx=0.5;
dz=0;
L=30;
im=L/dx+1;
hl=2;
hr=1;
out=60;
tpl=time/out;

%intitial condition
for i=1:im
    if i<im/2
        y(i)=hl;
        hf(i)=hl;
    else
        y(i)=hr;
        hf(i)=hr;
    end
end
   %VARIASI KONTRAKSI DAN EKSPANSI-------------
for i=1:im
    x(i)=(i-1)*dx;
    L1=8;  %dapat diganti *nb= L1+L2+L3=15
    L2=7;   %dapat diganti *nb= L1+L2+L3=15
    L3=15;  %dapat diganti *nb= L1+L2+L3=15
    B1=5;   %dapat diganti bebas
    B2=1;   %dapat diganti bebas
    B3=2;   %dapat diganti bebas
    if x(i)<=L1
        bmin=(B2-B1)/2;
        alfa=atand(bmin/L1);
        b(i)=x(i)*tand(alfa)*2+B1;
    elseif x(i)>L1 && x(i)<=(L1+L2)
        b(i)=B2;     
    else
        bmin2=(B2-B3)/2;
        alfa2=atand(bmin2/L3);
        b(i)=B2-(x(i)-L3)*tand(alfa2)*2;
    end
    if i<=11
        b(i)=1;
        n(i)=0.01;
    elseif i<=21
        b(i)=3;
        n(i)=0.01;
    elseif i<=31
        b(i)=


    z(i)=-x(i)*So;
    Q(i)=0;
    Qf(i)=0;
    A(i)=b(i)*y(i);
    Af(i)=b(i)*y(i);
    p(i)=b(i)/2+B1/2;
    q(i)=b(i)/2*-1+B1/2;

end

subplot(2,1,2);
plot(x,p,'red',LineWidth=3);
hold on
plot(x,q,'red',LineWidth=3);
title('Tampak Atas Saluran');
xlabel('Jarak (m)');
ylabel('Lebar (m)');
grid on
grid minor

% --------------------- E -- N -- G -- I -- N -- E-- ---------------------
for t=1:tm
%=================================FTCS====================================
    %kontinuitas++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    for i=2:im-1
        Ab(i)=Af(i)-dt/2/dx*(Qf(i+1)-Qf(i-1));
    end
    %momentum ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    for i=2:im-1
        Sof=-dz/dx;
        uf=Qf(i)/Af(i);
        hf(i)=Ab(i)/b(i);
        R=Af(i)/(b(i)+2*hf(i));
        Sf=n(i)^2*abs(uf)*uf/R^(4/3);
        bro1=dt/2/dx*(Qf(i+1)^2/Af(i+1)-Qf(i-1)^2/Af(i-1));
        bro2=g*Af(i)*dt*(Sof-Sf);
        bro3=g*Af(i)*dt/2/dx*(hf(i+1)-hf(i-1));
        Qb(i)=Qf(i)-bro1+bro2-bro3;
    end

    % syarat batas +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Ab(1)=Ab(2);
    Qb(1)=0;
    Ab(im)=Ab(im-1);
    Qb(im)=0;
    hf(1)=hf(2);
    hf(im)=hf(im-1);

    Af=Ab;
    Qf=Qb;

    %FILTER HANSEN +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Afill=Af;
    Qfill=Qf;
    Cofil=0.99;
    for i=2:im-1
        Af(i)=Cofil*Afill(i)+(1-Cofil)*(Afill(i-1)+Afill(i+1))/2;
        Qf(i)=Cofil*Qfill(i)+(1-Cofil)*(Qfill(i-1)+Qfill(i+1))/2;

    end

%==============================LAX WENDROFF================================


    %STEP 1----------------------------------------------------------------
    for i=1:im-1
        %kontinuitas
        A1=0.5*(A(i+1)+A(i)); % At, i+1/2
        Ahalf(i)=A1-dt/2/dx*(Q(i+1)-Q(i)); % At+1, i+1/2
        yhalf(i) =Ahalf(i)/b(i); %yt+1, i+1/2
        %Momentum
        So      =-dz/dx;
        Q1      =0.5*(Q(i)+Q(i+1)); % Qt, i+1/2;
        u1      =Q1/A1;
        R1      =A1/(b(i)+2*(A1/b(i)));
        Sf1     =n(i)^2*g*abs(u1)*u1/R1^(4/3);
        lax1    =dt/2/dx*(Q(i+1)^2/A(i+1)-Q(i)^2/A(i));
        lax2    =g*A1*dt/2*((y(i+1)-y(i))/dx+Sf1-So);
        Qhalf(i)=Q1-lax1-lax2; %Qt+1/2, i+1/2
    end

    %STEP 2----------------------------------------------------------------
    for i=2:im-1
        %kontinuitas
        A2(i)   =A(i)-dt/dx*(Qhalf(i)-Qhalf(i-1));    %A t+1, i
        y2(i)   =A2(i)/b(i);                          %y t+1, i
        %Momentum
        So      =-dz/dx;
        u2      =Q(i)/A(i);
        R2      =A(i)/(b(i)+2*y(i));
        Sf2     =n^2*g*abs(u2)*u2/R2^(4/3);
        wen1    =dt/dx*(Qhalf(i)^2/Ahalf(i)-Qhalf(i-1)^2/Ahalf(i-1));
        wen2    =g*A(i)*dt*((yhalf(i)-yhalf(i-1))/dx+Sf2-So);
        Qb(i)   =Q(i)-wen1-wen2;
    end

    %boundary condition----------------------------------------------------
    A2(1)=A2(2);
    Qb(1)=0;
    A2(im)=A2(im-1);
    Qb(im)=0;
    y2(1)=y2(2);
    y2(im)=y2(im-1);
    A=A2;
    Q=Qb;
    y=y2;

    %FILTER HANSEN---------------------------------------------------------
    Afil=A;
    Qfil=Q;
    Cofil=0.99;
    for i=2:im-1
        A(i)=Cofil*Afil(i)+(1-Cofil)*(Afil(i-1)+Afil(i)+Afil(i+1))/3;
        Q(i)=Cofil*Qfil(i)+(1-Cofil)*(Qfil(i-1)+Qfil(i)+Qfil(i+1))/3;
        u(i)=Q(i)/A(i);
        u=u1;
    end

    %=====================================================================
    %ANALITIK
    ts=t*dt;
    hm=(hl+hr)/2;
    cm=sqrt(hm*g);
    xa=(im-1)/2*dx-ts*(g*hl)^0.5;
    xb=(im-1)/2*dx+ts*(2*(g*hl)^0.5-3*cm);
    xc=(im-1)/2*dx+ts*2*cm^2*((g*hl)^0.5-cm)/(cm^2-g*hr);
    for i=1:im
        if (i-1)*dx<=xa
            han(i)=hl;
        elseif (i-1)*dx<=xb
            han(i)=4/9/g*((g*hl)^0.5-((i*dx-im/2*dx)/(2*ts)))^2;
        elseif (i-1)*dx<=xc
            han(i)=hm;
        elseif(i-1)*dx>=xc
            han(i)=hr;
        end
    end

    %PLOTTING AKTIVKAN SALAH SATU (GRAFIK ATAU ANIMASI)
%     plotting grafik
%         CN=(dt/dx);
%         if mod(t*dt,tplot)<=0
%             hold on
%             yyaxis left
%             txt4=['CN=',num2str(CN)];
%             txt1=['Numerik t= ',num2str(t*dt), ' jam'];
%             plot(x,(y+z),'blue','DisplayName',txt1);
%             txt2=['Analitik t= ',num2str(t*dt), ' jam'];
%             plot(x,han,'red','DisplayName',txt2);
%             ylim([0.9 1.6])
%             ylabel('Tinggi Muka Air (m)')
%             hold on
%             yyaxis right
%             ylim([-0.1 1.6])
%             xlabel('Jarak (m)');
%             title(['Hasil Simulasi Dam Break dengan Lax Wendroff, Cf=',num2str(Cofil)],txt4,FontSize=16)
%             legend("show" ,Location="bestoutside")
%             hold off
%         end


    %plotting ANIMASI
    if mod(t*dt,tpl)<=0
        subplot(2,1,1)
        p1=plot(x,y,'blue',LineWidth=1);
        hold on
        p2=plot(x,hf,'magenta');
        p3=plot(x,han,'black');
        axis([0 L 0.5 (hl+1)]);
        title(['LaxWendroff vs FTCS  t= ' num2str(ts)]);
        legend([p1,p2,p3],'Lax Wendroff','FTCS','Analitik',Location='bestoutside');
        xlabel('Jarak (m)')
        ylabel('Muka air (m)');
        pause(0.01);
        hold off;
    end

end

