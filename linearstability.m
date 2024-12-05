%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROGRAM linearstability.m
% THIS PROGRAM FINDS THE UNSTABLE WAVELENGHTS OF A MEAN VELOCITY PROFILE
% COUPLED TO THE MOTION OF CROPS FOR INCREMENTING REDUCED VELOCITIES.
% IT RETURNS THE CURVES OF THE MOST UNSTABLE WAVELENGTH VERSUS UR
%
% THIS PROGRAM PRODUCES THE PLOT OF FIGURE 
%
% MORE INFORMATION CAN BE FOUND IN THE PAPER "MODELLING WAVING CROPS USING
% LARGE EDDY SIMULATION" BY S. DUPONT, F. GOSSELIN, C. PY. E DE LANGRE, P.
% HEMON AND Y. BRUNET. SEE ALSO THE THESIS BE F. GOSSELIN "MECANISMES
% D'INTERACTIONS FLUIDE-STRUCTURE ENTRE ECOULEMENTS ET VEGETATION."
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linsty='-b'
linsty2='-g'
Nk=50;                          %RESOLUTION IN WAVELENGTH. NUMBER OF WAVELENGTH/WAVENUMBER
lammin=2;                       %MINIMAL WAVELENGTH (MEASURED IN h)
lammax=8;                       %MAXIMAL WAVELENGTH (MEASURED IN h)
dlam=(lammax-lammin)/(Nk-1);    %WAVELENGTH INCREMENT
Ny=400;                         %RESOLUTION IN THE VERTICAL DIRECTION. NUMBER OF POINTS IN THE FINITE DIFFERENCE
                                %SHOULD BE ABOVE 150 FOR A NUMERICALLY
                                %REASONABLE ANSWER AND NOT MORE THAN A
                                %1000 FOR A REASONABLE COMPUTATION TIME.
H=12;                           %HEIGHT OF THE DOMAIN (MEASURED IN h)
dy=(H-0)/Ny;
Rey=-1                          %TURBULENT REYNOLDS NUMBER. WHEN THE VALUE IS -1, THE TURBULENT PROFILE IS DEFINED IN THE FILE invreyfunc.m
hc=0.69;                        %CANOPY HEIGHT (m)
fc=1.05;                        %NATURAL FREQUENCY OF THE PLANTS (hz)
M=0.0046/1.23/hc/0.05^2;        %MASS NUMBER
xi=0.0879;                      %DAMPING COEFFICIENT



%REDUCED VELOCITY RANGE. NOTE THAT UR DEFINED HERE DIFFERS FROM THAT IN THE
%PAPER FROM A FACTOR 2pi
allU=[0.1 0.5:0.5:10]/2/pi



clear y Ub ReyYD CUYD
for n1=1:Ny
    y(n1,1)=n1*dy;
    Ub(n1,:)=Ubfunc(y(n1,1),H);
    invreyYD(n1,:)=invreyfunc(y(n1,1),Rey,H);
    C1(1,:)=cfunc(y(n1,1),H);
    CUYD(n1,1)=C1(1,1)*Ub(n1,1);                    %product of C times U
    CUYD(n1,2)=C1(1,2)*Ub(n1,1)+C1(1,1)*Ub(n1,2);   %derivative of (product of C times U)
end


timezero=clock;
clear maxvsre nuc
nuc=0;
for UR=allU %REDUCED VELOCITY LOOP
    
    Ucanopytop=UR*hc*fc*2*pi;
    
    nuc=nuc+1;
    
    clear maxrootvsk fourfreqs maxrootsolid1 fourfreqs
    for nk=1:Nk     %WAVENUMBER LOOP
        
        lam=lammin+dlam*(nk-1);
        k=2*pi/lam;
        
        valandvec=findomegacanopymanyfreqs(Ny,dy,Ub,k,invreyYD,CUYD,M,UR,xi);
        omk=valandvec(1,1)+i*valandvec(1,2);
        
        [nfreqs,nout]=size(valandvec);

        maxrootvsk(nk,1)=k;                             %WAVENUMBER
        maxrootvsk(nk,2)=real(omk);                     %REAL FREQUENCY
        maxrootvsk(nk,3)=imag(omk);                     %GROWTH RATE
        maxrootvsk(nk,4)=2*pi/k;                        %WAVELENGTH
        maxrootvsk(nk,5)=real(omk)/Ucanopytop*hc;     %DIMENSIONAL FREQUENCY OF THE INSTABILITY
        maxrootvsk(nk,6)=valandvec(1,nout);             %PERCENTAGE OF ALL ENERGY THAT IS SOLID
        
        
        %THIS NEXT PIECE OF CODE IS A REMNANT OF MY IDEA TO LOOK AT MODES
        %WITH MORE ENERGY IN THE SOLID. I LEFT IT THERE FOR MY POSSIBLE
        %THOUGH IMPROBABLE FUTURE REUSE.
        %{
        [garbage,nsol]=max(valandvec(2:3,2));
        maxrootsolid1(nk,1)=k;
        maxrootsolid1(nk,2)=valandvec(nsol+1,1);
        maxrootsolid1(nk,3)=valandvec(nsol+1,2);
        maxrootsolid1(nk,4)=2*pi/k;
        maxrootsolid1(nk,5)=valandvec(nsol+1,1)/Ucanopytop*hc;
                
        fourfreqs(nk,1)=k;
        fourfreqs(nk,2)=2*pi/k;
        fourfreqs(nk,3)=real(omk)/Ucanopytop*hc;
        for nftp=1:nfreqs
            fourfreqs(nk,1+3*nftp)=valandvec(nftp,1);
            fourfreqs(nk,2+3*nftp)=valandvec(nftp,2);
            fourfreqs(nk,3+3*nftp)=valandvec(nftp,nout);
        end
        %}
        
        
        percentdone=(nuc-1+(nk)/Nk)/(size(allU,2))*100
        timeelapsed=etime(clock,timezero);
        timeremaining=timeelapsed/percentdone*(100-percentdone)/60
        
    end
    
      
    [garbage,nmaxk]=max(maxrootvsk(:,3));
    temprey=invreyfunc(1,Rey,H);
    maxvsre(nuc,1)=1/temprey(1,1);                                      %REYNOLD NUMBER AT CANOPY TOP
    maxvsre(nuc,2)=maxrootvsk(nmaxk,1);                                 %k, WAVENUMBER
    maxvsre(nuc,3)=maxrootvsk(nmaxk,2);                                 %omega_r, REAL FREQUENCY
    maxvsre(nuc,4)=maxrootvsk(nmaxk,3);                                 %omega_i, GROWTH RATE
    maxvsre(nuc,5)=maxrootvsk(nmaxk,4);                                 %lambda/h, WAVELENGHT
    maxvsre(nuc,6)=maxrootvsk(nmaxk,5);                                 %DIMENSIONAL FREQUENCY OF THE INSTABILITY
    maxvsre(nuc,7)=Ucanopytop;                                          %U @ y=h
    maxvsre(nuc,8)=UR;                                                  %UR
    maxvsre(nuc,9)=UR*maxrootvsk(nmaxk,2);                              %UR*freq
    maxvsre(nuc,10)=maxrootvsk(nmaxk,6);                                %PERCENTAGE OF ALL ENERGY THAT IS SOLID
    maxvsre(nuc,11)=maxrootvsk(nmaxk,2)/maxrootvsk(nmaxk,1);            %Ucp/Uh PHASE VELOCITY

    
    %THIS NEXT PIECE OF CODE IS A REMNANT OF MY IDEA TO LOOK AT MODES
    %WITH MORE ENERGY IN THE SOLID. I LEFT IT THERE FOR MY POSSIBLE
    %THOUGH IMPROBABLE FUTURE REUSE.
    %{
    [garbage,nmaxk2]=max(maxrootsolid1(:,3));
    temprey=invreyfunc(1,Rey,H,Ucanopytop);
    maxvsre2(nuc,1)=1/temprey(1,1);                                      %Reynolds
    maxvsre2(nuc,2)=maxrootsolid1(nmaxk2,1);                                 %k
    maxvsre2(nuc,3)=maxrootsolid1(nmaxk2,2);                                 %omega_r
    maxvsre2(nuc,4)=maxrootsolid1(nmaxk2,3);                                 %omega_i
    maxvsre2(nuc,5)=maxrootsolid1(nmaxk2,4);                                 %lambda/h
    maxvsre2(nuc,6)=maxrootsolid1(nmaxk2,5);                                 %?
    maxvsre2(nuc,7)=Ucanopytop;                                          %U @ y=h
    maxvsre2(nuc,8)=UR;                                                  %UR
    maxvsre2(nuc,9)=Ucanopytop/hc/(fc*2*pi)*maxrootsolid1(nmaxk2,2);     %UR*freq
    %}
    
    %PLOTS THE UNSTABLE CURVE AT THE GIVEN REDUCED VELOCITY
    figure(3)
    subplot(3,1,1)
    hold on
    plot(maxrootvsk(:,4),maxrootvsk(:,2),linsty)
    %plot(maxrootsolid1(:,4),maxrootsolid1(:,2),linsty2)
    ylabel('\omega_r')
    subplot(3,1,2)
    hold on
    plot(maxrootvsk(:,4),maxrootvsk(:,3),linsty)
    %plot(maxrootsolid1(:,4),maxrootsolid1(:,3),linsty2)
    ylabel('\omega_i')
    subplot(3,1,3)
    hold on
    plot(maxrootvsk(:,4),maxrootvsk(:,6),linsty)
    %plot(maxrootsolid1(:,4),maxrootsolid1(:,3),linsty2)
    ylabel('percent E-solid')
    xlabel('\lambda/h')
   
    
    markerbank=[' -ok';' -xb';' -+r';'--xk'];
    
    %THIS NEXT PIECE OF CODE IS A REMNANT OF MY IDEA TO LOOK AT MODES
    %WITH MORE ENERGY IN THE SOLID. I LEFT IT THERE FOR MY POSSIBLE
    %THOUGH IMPROBABLE FUTURE REUSE.
    %{
    figure(3)
    for nftp=1:nfreqs
        subplot(3,1,1)
        hold on
        plot(fourfreqs(:,1),fourfreqs(:,1+nftp*3)*UR,markerbank(nftp,:))
        ylabel('\omega_r U_R')
        if nftp==1
            subplot(3,1,2)
            hold on
            plot(fourfreqs(:,1),fourfreqs(:,2+nftp*3),markerbank(nftp,:))
            ylabel('\omega_i')
        else
            subplot(3,1,3)
            hold on
            plot(fourfreqs(:,1),fourfreqs(:,2+nftp*3),markerbank(nftp,:))
            ylabel('\omega_i')
            xlabel('k')
        end
    end
   %}
    %{
    figure(nuc+1)
    for nftp=1:nfreqs
        subplot(4,1,1)
        hold on
        title(strcat('Uh=',num2str(Ucanopytop)))
        plot(fourfreqs(:,2),fourfreqs(:,1+nftp*3)*UR,markerbank(nftp,:))
        ylabel('\omega_r U_R')
        xlim([1 15])
        if nftp==1
            subplot(4,1,2)
            hold on
            plot(fourfreqs(:,2),fourfreqs(:,2+nftp*3),markerbank(nftp,:))
            ylabel('\omega_i')
            xlim([1 15])
        else
            subplot(4,1,3)
            hold on
            plot(fourfreqs(:,2),fourfreqs(:,2+nftp*3),markerbank(nftp,:))
            ylabel('\omega_i')
            xlim([1 15])

        end

        subplot(4,1,4)
        plot(fourfreqs(:,2),fourfreqs(:,3+nftp*3),markerbank(nftp,:))
        hold on
        ylabel('E-fraction')
        xlabel('\lambda /h')
        xlim([1 15])

    end
%}

    %SAVE THE DATA FOR A GIVEN REDUCED VELOCITY. TO USE IT, UNCOMMENT THE
    %CSVWRITE LINE
    str1='maxrootvsk';
    str2='-Uh=';
    str3=num2str(Ucanopytop);
    strend='.csv';
    filename=strcat(str1,str2,str3,strend);
    %csvwrite(filename,maxrootvsk,0,0);
    
    %THIS NEXT PIECE OF CODE IS A REMNANT OF MY IDEA TO LOOK AT MODES
    %WITH MORE ENERGY IN THE SOLID. I LEFT IT THERE FOR MY POSSIBLE
    %THOUGH IMPROBABLE FUTURE REUSE.
    str1='threefreqs';
    str2='-M=';
    str3=num2str(M);
    str4='-Ny=';
    str5=num2str(Ny);
    str6='-Nk=';
    str7=num2str(Nk);
    str8='-Re=';
    str9=num2str(Rey);
    str10='-UR=';
    str11=num2str(UR);
    str12='-xi=';
    str13=num2str(xi);
    strend='.csv';
    filename=strcat(str1,str2,str3,str4,str5,str6,str7,str8,str9,str10,str11,str12,str13,strend);
    %csvwrite(filename,fourfreqs,0,0);

    imomega1=interp1(maxrootvsk(:,1),maxrootvsk(:,3),1,'nearest');    
end


%THIS SAVES THE FREQ, k, lambda, ETC VERSUS UR CURVES
str1='final-modefluide-';
str2='-M=';
str3=num2str(M);
str4='-Ny=';
str5=num2str(Ny);
str6='-Nk=';
str7=num2str(Nk);
str8='-Re=';
str9=num2str(Rey);
str10='-URmax=';
str11=num2str(UR);
str12='-xi=';
str13=num2str(xi);
strend='.csv';
filename=strcat(str1,str2,str3,str4,str5,str6,str7,str8,str9,str10,str11,str12,str13,strend);
%csvwrite(filename,maxvsre,0,0);



%PRODUCES THE PLOT OF FIGURE 16
figure(16)
subplot(4,1,1)
hold on
plot(2*pi*maxvsre(:,8),maxvsre(:,5),'-k')
ylabel('lambdaoh')
xlim([0 10])
ylim([3 6])
box on
subplot(4,1,2)
hold on
clear pl lege
pl(1)=plot(2*pi*maxvsre(:,8),maxvsre(:,9),'-k');
xlim([0 10])
ylim([0 2])
box on
ylabel('fof0')
subplot(4,1,3)
hold on
plot(2*pi*maxvsre(:,8),maxvsre(:,4),'-k')
ylabel('omegai')
xlim([0 10])
ylim([0 0.1])
%xlabel('UR')
box on
subplot(4,1,4)
hold on
plot(2*pi*maxvsre(:,8),maxvsre(:,10),'-k')
ylabel('eta')
xlim([0 10])
%ylim([0 0.1])
xlabel('UR')
box on

%PRODUCES THE PLOT OF FIGURE 13
figure(13)
subplot(3,1,1)
hold on
plot(2*pi*maxvsre(:,8),maxvsre(:,5),'-k')
ylabel('lambdaoh')
xlim([0 10])
ylim([3 6])
box on
subplot(3,1,2)
hold on
clear pl lege
pl(1)=plot(2*pi*maxvsre(:,8),maxvsre(:,9),'-k');
xlim([0 10])
ylim([0 2])
box on
ylabel('fof0')
subplot(3,1,3)
hold on
plot(2*pi*maxvsre(:,8),maxvsre(:,11),'-k')
ylabel('U_{cp}/U_h')
xlim([0 10])
%ylim([0 0.1])
xlabel('UR')
box on