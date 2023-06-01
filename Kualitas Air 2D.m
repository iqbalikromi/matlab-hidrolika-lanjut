%UAS-PEMODELAN KUALITAS AIR - 2022 
%AUTHOR:: AKHMAD IQBAL IKROMI - 25021042

clear
clf
clc

time=360*2;
dt=0.2;
loop=time/dt;
% output=time/(3600/5);
output=time;
L=5500; %m
b=7000; %m
n=0.03;
g=9.81;


dx=100;
dy=100;

sox=0;
soy=0;

im=L/dx+1;
jm=b/dy+1;
ac=200;

Pe=2;
%initial condition
for i=1:im
    for j=1:jm
        x(i,j)=(i-1)*dx;
        y(i,j)=(j-1)*dy;
        u(i,j)=0;
        Ltot(i,j)=i;
        Btot(i,j)=j;
        v(i,j)=0;
        h(i,j)=3;
        c(i,j)=0;
        z(i,j)=0;
    end
end

%lokasi inlet dan outlet
%coding ini menggunakan dx=100 sehingga dalam penentuan nilai i -nya
%berbeda-beda tiap kondisi, yang ditulis sebagai berikut.

%INLET 1
in1awalj=1700/dy+1;
in1akhirj=1900/dy; 
in1i=0/dx+1;
for i=in1i:3
    for j=in1awalj:in1akhirj
        ilet1i(i,j)=(i-1)*dx;
        ilet1j(i,j)=(j-1)*dy;
    end
end


%INLET 2
in2awalj=3700/dy+1;
in2akhirj=3700/dy+1; 
in2i=0/dx+1;
for i=in2i:3
    for j=in2awalj:in2akhirj
        ilet2i(i,j)=(i-1)*dx;
        ilet2j(i,j)=(j-1)*dy;
    end
end

%INLET3
in3awali=4000/dx+1;
in3akhiri=4000/dx+1;
in3j=7000/dx+1;
for i=in3awali:in3akhiri
    for j=(in3j-3):in3j
        ilet3i(i,j)=(i-1)*dx;
        ilet3j(i,j)=(j-1)*dy;
    end
end

%OUTLET
outawalj=0/dy+1;
outakhirj=300/dy+1;
outi=5500/dx+1;
for i=outi-3:outi
    for j=outawalj:outakhirj
        iout(i,j)=(i-1)*dx;
        jout(i,j)=(j-1)*dy;
    end
end


%LOKASI DARATAN
%DARATAN 1 ukuran 900 x 2200
d1awali=0/dx+1;
d1akhiri=2200/dx+1;
d1awalj=0/dy+1;
d1akhirj=900/dy+1;
for i=d1awali:d1akhiri
    for j=d1awalj:d1akhirj
        B1i(i,j)=(i-1)*dx;
        B1j(i,j)=(j-1)*dy;
    end
end

%DARATAN 2 ukuran 3200 x 2000
d2awali=0/dx+1;
d2akhiri=2000/dx+1;
d2awalj=3800/dy+1;
d2akhirj=7000/dy+1;
for i=d2awali:d2akhiri
    for j=d2awalj:d2akhirj
        B2i(i,j)=(i-1)*dx;
        B2j(i,j)=(j-1)*dy;
    end
end

%DARATAN 3 ukuran 950 x 3200
d3awali=4500/dx+1;
d3akhiri=5500/dx+1;
d3awalj=3800/dy+1;
d3akhirj=7000/dy+1;
for i=d3awali:d3akhiri
    for j=d3awalj:d3akhirj
        B3i(i,j)=(i-1)*dx;
        B3j(i,j)=(j-1)*dy;
    end
end

%PLOTTING PRE-RUN 
% plot(x,y,'.','Color','blue');
% hold on
% plot(B1i,B1j,'.','Color','red');
% plot(B2i,B2j,'.','Color','red');
% plot(B3i,B3j,'.','Color','red');
% plot(ilet1i,ilet1j,'.','color','green');
% plot(ilet2i,ilet2j,'.','color','green');
% plot(ilet3i,ilet3j,'.','color','green');
% plot(iout,jout,'.','color','green');
% 
% axis equal
% title('Pre-Running Program Kualitas Air');
% pause
% clf

%UAS-PEMODELAN KUALITAS AIR - 2022 
%AUTHOR:: AKHMAD IQBAL IKROMI - 25021042

%engine
for t=1:loop
    %% HYDRODYNAMIC ENGINE 
    for i=1:im-1
        for j=1:jm-1
            %predictor

            %kontinuitas            
            udhdxp=u(i,j)*(h(i+1,j)-h(i,j))/dx;
            hdudxp=h(i,j)*(u(i+1,j)-u(i,j))/dx;
            vdhdxp=v(i,j)*(h(i,j+1)-h(i,j))/dy;
            hdvdxp=h(i,j)*(v(i,j+1)-v(i,j))/dy;
            hp(i,j)=h(i,j)-dt*(udhdxp+hdudxp+vdhdxp+hdvdxp);

            %momentum arah x
            sfx=n^2*u(i,j)*(u(i,j)^2+v(i,j)^2)^0.5/(h(i,j)^(4/3));
            mo1xp=u(i,j)*(u(i+1,j)-u(i,j))/dx;
            mo2xp=v(i,j)*(u(i,j+1)-u(i,j))/dy;
            mo3xp=g*(h(i+1,j)+z(i+1,j)-h(i,j)-z(i,j))/dx;
            mo4xp=g*(sox-sfx);
            up(i,j)=u(i,j)-dt*(mo1xp+mo2xp+mo3xp-mo4xp);

            %momentum arah y
            sfy=n^2*v(i,j)*(u(i,j)^2+v(i,j)^2)^0.5/(h(i,j)^(4/3));
            mo1yp=u(i,j)*(v(i+1,j)-v(i,j))/dx;
            mo2yp=v(i,j)*(v(i,j+1)-v(i,j))/dy;
            mo3yp=g*(h(i,j+1)+z(i,j+1)-h(i,j)-z(i,j))/dy;
            mo4yp=g*(soy-sfy);
            vp(i,j)=v(i,j)-dt*(mo1yp+mo2yp+mo3yp-mo4yp);     
        end
    end

    for i=2:im-1
        for j=2:jm-1

            %corrector
            uhalf=0.5*(u(i,j)+up(i,j));
            vhalf=0.5*(v(i,j)+vp(i,j));
            hhalf=0.5*(h(i,j)+hp(i,j));

            %kontinuitas
            udhdxc=uhalf*(hp(i,j)-hp(i-1,j))/dx;
            hdudxc=hhalf*(up(i,j)-up(i-1,j))/dx;
            vdhdxc=vhalf*(hp(i,j)-hp(i,j-1))/dy;
            hdvdxc=hhalf*(vp(i,j)-vp(i,j-1))/dy;
            hn(i,j)=hhalf-dt/2*(udhdxc+hdudxc+vdhdxc+hdvdxc);

            %momentum arah x
            sfx=n^2*uhalf*(uhalf^2+vhalf^2)^0.5/(hhalf^(4/3));
            mo1xc=uhalf*(up(i,j)-up(i-1,j))/dx;
            mo2xc=vhalf*(up(i,j)-up(i,j-1))/dy;
            mo3xc=g*(hp(i,j)+z(i,j)-hp(i-1,j)-z(i-1,j))/dx;
            mo4xc=g*(sox-sfx);
            un(i,j)=uhalf-dt/2*(mo1xc+mo2xc+mo3xc-mo4xc);

            %momentum arah y
            sfy=n^2*vhalf*(uhalf^2+vhalf^2)^0.5/(hhalf^(4/3));
            mo1yc=uhalf*(vp(i,j)-vp(i-1,j))/dx;
            mo2yc=vhalf*(vp(i,j)-vp(i,j-1))/dy;
            mo3yc=g*(hp(i,j)+z(i,j)-hp(i,j-1)-z(i,j-1))/dy;
            mo4yc=g*(soy-sfy);
            vn(i,j)=vhalf-dt/2*(mo1yc+mo2yc+mo3yc-mo4yc);

        end
    end

    %boundary condition

    %dinding kanan kiri
    for i=2:im-1
        %kanan
        hn(i,1)=hn(i,2);
        un(i,1)=un(i,2);
        vn(i,1)=0;
      
        %kiri
        hn(i,jm)=hn(i,jm-1);
        un(i,jm)=un(i,jm-1);
        vn(i,jm)=0;
    end
    %ujung kanan kiri
    for j=2:jm-1
        %tembok
        un(1,j)=0;
        un(im,j)=0;
        %ujung kiri
        hn(1,j)=hn(2,j);
        vn(1,j)=0;

        %ujung kanan
        hn(im,j)=hn(im-1,j);
        vn(im,j)=0;
    end

    %% SYARAT BATAS DARATAN
    %Daratan 1
    for i=d1awali:d1akhiri
        for j=d1awalj:d1akhirj
            un(i,j)=0;
            vn(i,j)=0;
        end
    end

    %Daratan 2
    for i=d2awali:d2akhiri
        for j=d2awalj:d2akhirj
            un(i,j)=0;
            vn(i,j)=0;
        end
    end

    %Daratan 3
    for i=d3awali:d3akhiri
        for j=d3awalj:d3akhirj
            un(i,j)=0;
            vn(i,j)=0;
        end
    end


    %% SYARAT BATAS SPILL POLUTAN DAN OUTLET MOMENTUM
    %inlet arah x
    %INLET 1
    for j=in1awalj:in1akhirj
        un(1,j)=2;
        hn(1,j)=3;
    end

    %INLET 2
    for j=in2awalj:in2akhirj
        un(1,j)=1;
        hn(1,j)=3;
    end

    %inlet arah y
    %INLET 3
    for i=in3awali:in3akhiri
        vn(i,jm)=-1.5;
        hn(i,jm)=3;
    end

    %OUTLET
    for j=outawalj:outakhirj
        vn(im,j)=0;
        un(im,j)=un(im-1,j); %BC Open
        hn(im,j)=3;
    end


    %bagian sudut

    hn(1,1)=(hn(1,2)+hn(2,1))/2;
    hn(1,jm)=(hn(2,jm)+hn(1,jm-1))/2;
    hn(im,1)=(hn(im-1,1)+hn(im,2))/2;
    hn(im,jm)=(hn(im,jm-1)+hn(im-1,jm))/2;

    un(1,1)=(un(1,2)+un(2,1))/2;
    un(1,jm)=(un(2,jm)+un(1,jm-1))/2;
    un(im,1)=(un(im-1,1)+un(im,2))/2;
    un(im,jm)=(un(im,jm-1)+un(im-1,jm))/2;

    vn(1,1)=(vn(1,2)+vn(2,1))/2;
    vn(1,jm)=(vn(2,jm)+vn(1,jm-1))/2;
    vn(im,1)=(vn(im-1,1)+vn(im,2))/2;
    vn(im,jm)=(vn(im,jm-1)+vn(im-1,jm))/2;

    %loop
    h=hn;
    v=vn;
    u=un;

    %filter
    hfl=h;
    vfl=v;
    ufl=u;
    cof=0.99;
    for i=2:im-1
        for j=2:jm-1
            h(i,j)=cof*hfl(i,j)+(1-cof)*(hfl(i+1,j)+hfl(i,j+1)+hfl(i,j)+hfl(i-1,j)+hfl(i,j-1))/5;
            v(i,j)=cof*vfl(i,j)+(1-cof)*(vfl(i+1,j)+vfl(i,j+1)+vfl(i,j)+vfl(i-1,j)+vfl(i,j-1))/5;
            u(i,j)=cof*ufl(i,j)+(1-cof)*(ufl(i+1,j)+ufl(i,j+1)+ufl(i,j)+ufl(i-1,j)+ufl(i,j-1))/5;
        end
    end

%UAS-PEMODELAN KUALITAS AIR - 2022 
%AUTHOR:: AKHMAD IQBAL IKROMI - 25021042
    %% WATER QUALITY ENGINE
    for i=1:im-2
        for j=1:jm-2
            
            %predictor

            suku1p=(u(i+1,j)*h(i+1,j)*c(i+1,j)-u(i,j)*h(i,j)*c(i,j))/dx;
            suku2p=(v(i,j+1)*h(i,j+1)*c(i,j+1)-v(i,j)*h(i,j)*c(i,j))/dy;
            suku3p=ac*h(i,j)*(c(i+2,j)-2*c(i+1,j)+c(i,j))/(2*dx)^2;
            suku4p=ac*h(i,j)*(c(i,j+2)-2*c(i,j+1)+c(i,j))/(2*dy)^2;
            cp(i,j)=(h(i,j)*c(i,j)+dt*(suku3p+suku4p-suku1p-suku2p))/h(i,j);
        end
    end
    for i=3:im-2
        for j=3:jm-2
            %corrector-final
            chalf=0.5*(c(i,j)+cp(i,j));
            suku1c=(u(i,j)*h(i,j)*c(i,j)-u(i-1,j)*h(i-1,j)*c(i-1,j))/dx;
            suku2c=(v(i,j)*h(i,j)*c(i,j)-v(i,j-1)*h(i,j-1)*c(i,j-1))/dy;
            suku3c=ac*h(i,j)*(c(i,j)-2*c(i-1,j)+c(i-2,j))/(2*dx)^2;
            suku4c=ac*h(i,j)*(c(i,j)-2*c(i,j-1)+c(i,j-2))/(2*dy)^2;
            cn(i,j)=(h(i,j)*chalf+dt*(suku3c+suku4c-suku1c-suku2c))/h(i,j);
        end
    end

    %syarat batas
        %bagian ujung kiri ketika j=3
    for j=3:jm-2
        cn(2,j)=cn(3,j);
        %bagian ujung kanan ketika i=im-1
        cn(im-1,j)=cn(im-2,j);
    end
        %dinding kanan i=3
    for i=3:im-2
        cn(i,2)=cn(i,3);
        %dinding kiri
        cn(i,jm-1)=cn(i,jm-2);
    end
        %bagian ujung kiri ketika i=2
    for j=2:jm-1
        cn(1,j)=cn(2,j);
        %bagian ujung kanan ketika i=im
        cn(im,j)=cn(im-1,j);
    end

        %dinding kanan ketiks j=2
    for i=2:im-1
        cn(i,1)=cn(i,2);
        %dinding kiri
        cn(i,jm)=cn(i,jm-1);
    end
    %% SYARAT BATAS INLET POLUTAN DAN OUTLET
    if (t*dt/3600)<=48
        %inlet arah x
        %INLET 1
        for j=in1awalj:in1akhirj
            cn(1,j)=10;
        end

        %INLET 2
        for j=in2awalj:in2akhirj
            cn(1,j)=5;
        end

        %inlet arah y
        %INLET 3
        for i=in3awali:in3akhiri
            cn(i,jm)=10;
        end
    else
              %inlet arah x
        %INLET 1
        for j=in1awalj:in1akhirj
            cn(1,j)=0;
        end

        %INLET 2
        for j=in2awalj:in2akhirj
            cn(1,j)=0;
        end

        %inlet arah y
        %INLET 3
        for i=in3awali:in3akhiri
            cn(i,jm)=0;
        end
    end


    


    %bagian sudut
    cn(1,1)=(cn(1,2)+cn(2,1))/2;
    cn(1,jm)=(cn(2,jm)+cn(1,jm-1))/2;
    cn(im,1)=(cn(im-1,1)+cn(im,2))/2;
    cn(im,jm)=(cn(im,jm-1)+cn(im-1,jm))/2;


    %Sisi Daratan 1
    for i=d1awali:d1akhiri
        for j=d1akhirj
            cn(i,j)=cn(i,j+1);            
        end
    end
    for i=d1akhiri
        for j=d1awalj:d1akhirj
            cn(i,j)=cn(i+1,j);
        end
    end


    %Daratan 2
    for i=d2awali:d2akhiri
        for j=d2awalj
            cn(i,j)=cn(i,j-1);    
        end
    end
        for i=d2akhiri
        for j=d2awalj:d2akhirj
            cn(i,j)=cn(i+1,j);    
        end
    end

    %Daratan 3
    for i=d3awali:d3akhiri
        for j=d3awalj
            cn(i,j)=cn(i,j-1);            
        end
    end
    for i=d3awali
        for j=d3awalj:d3akhirj
            cn(i,j)=cn(i-1,j);            
        end
    end

    % untuk bagian dalam daratan
    %Daratan 1
    for i=d1awali:d1akhiri-1
        for j=d1awalj:d1akhirj-1
            cn(i,j)=0;
        end
    end

    %Daratan 2
    for i=d2awali:d2akhiri-1
        for j=d2awalj+1:d2akhirj
            cn(i,j)=0;
        end
    end

    %Daratan 3
    for i=d3awali+1:d3akhiri
        for j=d3awalj+1:d3akhirj
            cn(i,j)=0;           
        end
    end

    %filter
    cfl=cn;
    cof=0.99;
    for i=2:im-1
        for j=2:jm-1
            cn(i,j)=cof*cfl(i,j)+(1-cof)*(cfl(i+1,j)+cfl(i,j+1)+cfl(i,j)+cfl(i-1,j)+cfl(i,j-1))/5;
        end
    end
    c=cn;

%UAS-PEMODELAN KUALITAS AIR - 2022 
%AUTHOR:: AKHMAD IQBAL IKROMI - 25021042
    %% plotting
    if mod(t*dt,output)==0
        figure(1)
        contour(x,y,c,100,ShowText="on");
        
        zticklabels({0 2 4 6 8 10})
   
        colormap("hsv")
     
        hold on
        quiver(x,y,u,v,5,'blue',ShowArrowHead="on",AutoScale="on");
        k1=plot(B1i,B1j,'x','Color','black');
        k2=plot(B2i,B2j,'x','Color','black');
        k3=plot(B3i,B3j,'x','Color','black');

        axis equal      

        xlim([0 L]);
        ylim([0 b]);
        xlabel('Panjang Danau (m)');
        ylabel('Lebar Danau (m)');

        title(['Penyebaran Polutan di Danau Saat t= ',num2str(t*dt/3600),' jam'])
        hold off
        saveas(gcf, ['t= ',num2str(t*dt/3600), 'jam.png']);

    end

end
%UAS-PEMODELAN KUALITAS AIR - 2022 
% AKHMAD IQBAL IKROMI - 25021042

