
% AUTHOR: AKHMAD IQBAL IKROMI
% COPY-PASTE NOT ALLOWED, READ AND LEARN ONLY

clear
clf
clc

simulasi= 12; %jam

% yn=5; %menggunakan data no 1

g=9.81;
m=0;
n=0.03; 
hlimit=0.5;
D=10;
% periode_pasut= 4000; %detik, simulasi minimal 3kali periode

datapasut=importdata("Pasut1.txt");
pasut=datapasut(:,2);
pasutx=transpose(pasut);
time=simulasi*3600;
dt=0.1;
tm=time/dt;
tpl=time/(simulasi*30);


infosal=importdata("PROGRAM\Gemoetry_KaliGung_OK.txt");

%data saluran

dxa=infosal(:,1);
dx=transpose(dxa);
L=sum(dx);
im=numel(dx);
for i=2:im
    x(1)=0;
    x(i)=x(i-1)+dx(i);
end

%lebar dasar
ba=infosal(:,4);
ma=infosal(:,5);
b=transpose(ba);
m=transpose(ma);


%dasar saluran
za=infosal(:,3);
eva=infosal(:,2);

z=transpose(za);
ev_tebing=transpose(eva);


for i=2:im
    dz(1)=0;
    dz(i)=z(i-1)-z(i);
    So(i)=dz(i)/dx(i);
    tsal(1)=ev_tebing(1)-z(1);
    tsal(i)=ev_tebing(i)-z(i);
end


%data hidrfograf
debit=importdata("PROGRAM\Debit_Q2.txt"); %pilih Q2 atau Q25
Qin=debit(:,2);
waktu=debit(:,1);

%initial condition

% yn=pasutx(1);
yn=2;

for i=1:im
     
    Q(i)=0;
    A(i)=(b(i)+m(i)*yn)*yn;
    h(i)=yn;
    p(i)=b(i)/2+b(1)/2;
    q(i)=b(i)/2*-1+b(1)/2;
    
end


%gambar saluran
for i=1:im
    for j=1:4
        %stasiun
        sta_Cross(i,1)=0;
        sta_Cross(i,2)=m(i)*tsal(i);
        sta_Cross(i,3)=sta_Cross(i,2)+b(i);
        sta_Cross(i,4)=sta_Cross(i,3)+sta_Cross(i,2);
        %elevasi
        ev_sal(i,1)=ev_tebing(i);
        ev_sal(i,2)=z(i);
        ev_sal(i,3)=z(i);
        ev_sal(i,4)=ev_tebing(i);
    end
    emin(i)=z(i)-0.3;
    if h(i)<ev_tebing(i)
        emax(i)=ev_tebing(i)+0.5;
    else
        emax(i)=h(i)+z(i)+0.5;
    end
    stamin(i)=0;
    stamax(i)=sta_Cross(i,4);

end

cawal=100;
cpol=300;
for i=1:im
    cf(i)=cawal;

    if i==10

        cf(i)=cpol;
    end
end
akhir=im+1;

for i=1:im
    for j=1:4
        hplot(i,j)=h(i)+z(i);
    end
end
for i=1:im
%     set(gcf, 'Position', get(0, 'Screensize'));
    subplot(9,5,i)

    area(sta_Cross(i,:),ev_sal(i,:),BaseValue=-5);
    axis([stamin(i) stamax(i) emin(i) emax(i)]);
    xlabel("Jarak (m)")

    hold on
    plot(sta_Cross(i,:),hplot(i,:));
    title(['Cross Section Sta. ',num2str(x(i)), ', t= 0 jam'])
end
subplot(9,5,akhir);
text(0.5,0.5,'PRESS ANY KEY TO RUN!', 'FontSize',18,Color='r')
axis off
pause
clf
%% program
for t=1:tm
    ts(t)=t*dt;
        %% Mac Cormack

%PREDICTOR
    for i=2:im-1
    %kontinuitas
        Ap(i)=A(i)-dt/dx(i)*(Q(i+1)-Q(i));
        hp1=(-b(i)-sqrt(b(i)^2+4*m(i)*Ap(i)))/(2*m(i));
        hp2=(-b(i)+sqrt(b(i)^2+4*m(i)*Ap(i)))/(2*m(i));
        if hp1>hp2
            hp(i)=hp1;
        else
            hp(i)=hp2;
        end

        if  hp(i)<=hlimit
            hp(i)=hlimit;
        else
            %Momentum
            u1      =Q(i)/A(i);
            Pp      =b(i)+2*h(i)*sqrt(m(i)^2+1);
            R1      =A(i)/Pp;
            Sf1     =n^2*g*abs(u1)*u1/R1^(4/3);
            cP1      =dt/dx(i)*(Q(i+1)^2/A(i+1)-Q(i)^2/A(i));
            cP2      =g*A(i)*dt*((h(i+1)-h(i))/dx(i)+Sf1-So(i));
            Qp(i)   =Q(i)-cP1-cP2;
        end
    end

    %boundary
    Ap(1)=Ap(2);
    hp(1)=hp(2);
    Qp(1)=Qp(2); 

%CORRECTOR
    for i=2:im-1
    %kontinuitas
        Ac(i)   =0.5*(Ap(i)+A(i))-dt/2/dx(i)*(Q(i)-Q(i-1));
        
        hc1=(-b(i)-sqrt(b(i)^2+4*m(i)*Ac(i)))/(2*m(i));
        hc2=(-b(i)+sqrt(b(i)^2+4*m(i)*Ac(i)))/(2*m(i));
        if hc1>hc2
            hc(i)=hc1;
        else
            hc(i)=hc2;
        end

        if  hc(i)<=hlimit
            hc(i)=hlimit;
        else
            %Momentum
            Ahalf   =0.5*(A(i)+Ap(i));
            hhalf   =0.5*(h(i)+hp(i));
            Qhalf   =0.5*(Q(i)+Qp(i));
            u2      =Qhalf/Ahalf;
            Pc      =b(i)+2*hhalf*sqrt(m(i)^2+1);
            R2      =Ahalf/Pc;
            Sf2     =n^2*g*abs(u2)*u2/R2^(4/3);
            C1      =dt/dx(i)/2*(Q(i)^2/A(i)-Q(i-1)^2/A(i-1));
            C2      =g*Ahalf*dt/2*((h(i)-h(i-1))/dx(i)+Sf2-So(i));
            Qb(i)   =Qhalf-C1-C2;
        end
    end
    %% WATER QUALITY
    %FTCS
    for i=2:im-1
        uf=Q(i)/A(i);
        adveksif(i)=dt*Q(i)/A(i)/2/dx(i)*(cf(i+1)-cf(i-1));
        Disnumf=D+dt/2*uf^2; %menggunakan dispersi numerik sesuai perubahan kecepatan (du)
        dispersif(i)=dt*(Disnumf*(cf(i+1)-2*cf(i)+cf(i-1)))/dx(i)/dx(i);
        cbf(i)=cf(i)-adveksif(i)+dispersif(i);
    end
    

%%
    %BC hulu
    cbf(1)=cf(1);
    cbf(im)=cbf(im-1);

    hc(1)=hc(2);
    Ac(1)=Ac(2);
    tintp=t*dt/3600;
    Qb(1)=interp1(waktu,Qin,tintp);

    %BC hilir
    Qb(im)=Qb(im-1);
%     hc(im)=yn+sin(1.5*pi*ts(t)/(3600*2));
    hc(im)=2+interp1(waktu,pasut,tintp);
    Ac(im)=(b(im)+m(im)*hc(im))*hc(im);

    A=Ac;
    Q=Qb;
    h=hc;
    u=Q/A;
    cf=cbf;


    %FILTER HANSEN untuk hidrodinamik
    Qfil=Q;
    hfil=h;
    cfil=cf;
    Cofil=0.995;
    for i=2:im-1

        h(i)=(Cofil*(hfil(i)+z(i))+(1-Cofil)*(hfil(i-1)+z(i-1)+hfil(i+1)+z(i+1))/2)-z(i);
        A(i)=(b(i)+m(i)*h(i))*h(i);
        Q(i)=Cofil*Qfil(i)+(1-Cofil)*(Qfil(i-1)+Qfil(i+1))/2;
        cf(i)=Cofil*cfil(i)+(1-Cofil)*(cfil(i-1)+cfil(i+1))/2;

    end


    wkt=round(tintp,2,"decimals");
    %plotting ANIMASI 
    if mod(t*dt,tpl)<=0
        subplot(2,1,1)
        yyaxis left
        figure(1);
        p1=plot(x,h+z,'blue',LineWidth=1);
        
        hold on
        p2=plot(x,z,'black--',LineWidth=0.9);
        
       
        axis([0 x(im) -5 20]);
        title({'Simulasi Penyebaran Polutan di Sungai Gung, Tegal'; ['t= ', num2str(wkt), ' Jam']});
        subtitle(['Konsentrasi Awal  Sungai= ', num2str(cawal),' ppm, dan Konsentrasi Polutan= ', num2str(cpol), ' ppm.']);
        
        xlabel('Jarak dari Hulu (m)')
        ylabel('Elevasi Muka Air (m)');

        yyaxis right
        hold off
        p3=plot(x,cf,'green',LineWidth=3,LineStyle=':');
        ylabel ('Konsentrasi Polutan (ppm)')
        legend([p1,p2,p3],'Muka Air','Bed Elv','Penyebaran Polutan');        
        ylim([100 120]);
        grid on
        grid minor

       
        subplot(2,2,3)
        
        w1=plot(waktu,Qin,'cyan',LineWidth=1.5);
        w2=xline(wkt,'red',LineWidth=2);
        title('Hidrograf Banjir Inflow Sungai Gung')
        ylabel('Debit (m3/dt)');
        xlabel('Waktu (jam)');
        legend([w1,w2],'Hidrograf','Waktu Simulasi');
        xlim([0 24]);
        grid on
        grid minor

        subplot(2,2,4)
        k1=plot(waktu,pasut);
        k2=xline(wkt,'red',LineWidth=2);
        title('Grafik Pasang Surut Muara Sungai Gung')
        ylabel('El. Muka Air (m)');
        xlabel('Waktu (jam)');
        xlim([0 24]);
        legend([k1,k2],'El. M.A','Waktu Simulasi',Location='best');
        grid on
        grid minor
        pause(0.05);
    end
end

pause
clf
%gambar saluran
for i=1:im

    emin(i)=z(i)-0.3;
    if h(i)<ev_tebing(i)
        emax(i)=ev_tebing(i)+0.5;
    else
        emax(i)=h(i)+z(i)+0.5;
    end
    stamin(i)=0;
    stamax(i)=sta_Cross(i,4);

end
for i=1:im
    for j=1:4
        hplot(i,j)=h(i)+z(i);
    end
end

for i=1:im

    subplot(9,5,i)

    area(sta_Cross(i,:),ev_sal(i,:),BaseValue=-5);
    axis([stamin(i) stamax(i) emin(i) emax(i)]);
    xlabel("Jarak (m)")

    hold on
    plot(sta_Cross(i,:),hplot(i,:));
    title(['Cross Section Sta. ',num2str(x(i)), ', t= ',num2str(wkt),'jam'])
end







