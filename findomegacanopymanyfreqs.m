
function answer=findomegacanopymanyfreqs(Ny,dy,Ub,k,invreyYD,CUYD,M,UR,xi)

%Damping coefficient
%xi=0.0879;%0.55;
%([A]-omega*[B]).v=0
%([B]^(-1).[A]-omega*[I]).v=0
clear A B
A=sparse(Ny,Ny);
B=sparse(Ny,Ny);


clear term1
for n2=1:Ny
    y1=n2*dy;
    CU=CUYD(n2,1);
    dCU=CUYD(n2,2);
    
    if y1<2
        chi=y1;
        dchi=1;
        term1(n2,1)=(chi^2*CU);
        B(2,n2+2)=-dy/M/k*(CU*dchi+dCU*chi);
        NN2=n2;
    end
end
A(1,1)=1;       %q point point
B(1,2)=1;       %q point
A(2,1)=-i*2*xi/UR;
A(2,2)=UR^(-2);
B(2,1)=1;
B(2,2)=sum(term1,'double')*i*dy/M;

for n1=1:Ny
    y1=n1*dy;

    if y1<2
        chi=y1;
        dchi=1;
    else
        chi=0;
        dchi=0;
    end
    
    invrey=invreyYD(n1,1);
    invdrey=invreyYD(n1,2);
    invddrey=invreyYD(n1,3);
    
    CU=CUYD(n1,1);
    dCU=CUYD(n1,2);
    
    if n1==1
%{
        A(1,1+4)=(-1*invrey)*(-1)/dy^4;
        A(1,1+3)=(-1*invrey)*(6)/dy^4;
        A(1,1+2)=(-1*invrey)*(-14)/dy^4;
        A(1,1+1)=(-1*invrey)*(16)/dy^4+(i*k*Ub(n1,1)+2*k^2*invrey)/dy^2;
        A(1,1+0)=(-1*invrey)*(-9)/dy^4-(i*k*Ub(n1,1)+2*k^2*invrey)*2/dy^2-i*k^3*Ub(n1,1)-i*k*Ub(n1,3)-k^4*invrey;
        
        B(1,1+1)=i/dy^2;
        B(1,1+0)=-2*i/dy^2-i*k^2;
%}
        %{
        %VERSION 1 WITH PSI''''=0 @y=0
        A(n1+2,2)=-(dCU*chi+CU*dchi);
        A(n1+2,n1+2+2)=-1/dy^4*invrey;
        A(n1+2,n1+2+1)=4/dy^4*invrey+(i*k*Ub(n1,1)+2*k^2*invrey)/dy^2;
        A(n1+2,n1+2+0)=-6/dy^4*invrey-(i*k*Ub(n1,1)+2*k^2*invrey)*2/dy^2-i*k^3*Ub(n1,1)-i*k*Ub(n1,3)-k^4*invrey+1/dy^4*invrey;
        %A(n1,n1-1)=4/dy^4*invrey+(i*k*Ub(n1,1)+2*k^2*invrey)/dy^2;
        %A(n1,n1-2)=-1/dy^4*invrey;
        %}
        %VERSION 2
        
        A(n1+2,2)=i*k*(dCU*chi+CU*dchi);
        A(n1+2,n1+2+3)=(-1/10)*(-1/dy^4*invrey+2/2/dy^3*invdrey);
        A(n1+2,n1+2+2)=-1/dy^4*invrey-2/2/dy^3*invdrey+(3/5)*(-1/dy^4*invrey+2/2/dy^3*invdrey);
        A(n1+2,n1+2+1)=4/dy^4*invrey+4/2/dy^3*invdrey+(i*k*Ub(n1,1)+2*k^2*invrey+1*invddrey+CU)/dy^2+(4*k^2*invdrey+dCU)/2/dy+(-7/5)*(-1/dy^4*invrey+2/2/dy^3*invdrey);
        A(n1+2,n1+2+0)=-6/dy^4*invrey-(i*k*Ub(n1,1)+2*k^2*invrey+1*invddrey+CU)*2/dy^2-i*k^3*Ub(n1,1)-i*k*Ub(n1,3)-k^4*invrey+k^2*invddrey+(2/5)*(-1/dy^4*invrey+2/2/dy^3*invdrey);
        %A(n1,n1-1)=4/dy^4*invrey+(i*k*Ub(n1,1)+2*k^2*invrey)/dy^2;
        %A(n1,n1-2)=-1/dy^4*invrey;
        
        B(n1+2,n1+2+1)=i/dy^2;
        B(n1+2,n1+2+0)=-2*i/dy^2-i*k^2;
        %B(n1,n1-1)=i/dy^2;
        
    elseif n1==2
        A(n1+2,2)=i*k*(dCU*chi+CU*dchi);
        A(n1+2,n1+2+2)=-1/dy^4*invrey-2/2/dy^3*invdrey;
        A(n1+2,n1+2+1)=4/dy^4*invrey+4/2/dy^3*invdrey+(i*k*Ub(n1,1)+2*k^2*invrey+1*invddrey+CU)/dy^2+(4*k^2*invdrey+dCU)/2/dy;
        A(n1+2,n1+2+0)=-6/dy^4*invrey-(i*k*Ub(n1,1)+2*k^2*invrey+1*invddrey+CU)*2/dy^2-i*k^3*Ub(n1,1)-i*k*Ub(n1,3)-k^4*invrey+k^2*invddrey;
        A(n1+2,n1+2-1)=4/dy^4*invrey-4/2/dy^3*invdrey+(i*k*Ub(n1,1)+2*k^2*invrey+1*invddrey+CU)/dy^2-(4*k^2*invdrey+dCU)/2/dy;
        %A(n1,n1-2)=-1/dy^4*invrey+2/2/dy^3*invdrey;
        
        B(n1+2,n1+2+1)=i/dy^2;
        B(n1+2,n1+2+0)=-2*i/dy^2-i*k^2;
        B(n1+2,n1+2-1)=i/dy^2;
    elseif n1==Ny
        %{
        %VERSION 1 WITH PSI''''=0 @y=0
        %A(n1,n1+2)=-1/dy^4*invrey;
        %A(n1,n1+1)=4/dy^4*invrey+(i*k*Ub(n1,1)+2*k^2*invrey)/dy^2;
        A(n1+2,n1+2+0)=-6/dy^4*invrey-(i*k*Ub(n1,1)+2*k^2*invrey)*2/dy^2-i*k^3*Ub(n1,1)-i*k*Ub(n1,3)-k^4*invrey+1/dy^4*invrey;
        A(n1+2,n1+2-1)=4/dy^4*invrey+(i*k*Ub(n1,1)+2*k^2*invrey)/dy^2;
        A(n1+2,n1+2-2)=-1/dy^4*invrey;
        %}
        %VERSION 2
        %A(n1,n1+2)=-1/dy^4*invrey;
        %A(n1,n1+1)=4/dy^4*invrey+(i*k*Ub(n1,1)+2*k^2*invrey)/dy^2;
        A(n1+2,n1+2+0)=-6/dy^4*invrey-(i*k*Ub(n1,1)+2*k^2*invrey+1*invddrey+CU)*2/dy^2-i*k^3*Ub(n1,1)-i*k*Ub(n1,3)-k^4*invrey+k^2*invddrey+(2/5)*(-1/dy^4*invrey-2/2/dy^3*invdrey);
        A(n1+2,n1+2-1)=4/dy^4*invrey-4/2/dy^3*invdrey+(i*k*Ub(n1,1)+2*k^2*invrey+1*invddrey+CU)/dy^2-(4*k^2*invdrey+dCU)/2/dy+(-7/5)*(-1/dy^4*invrey-2/2/dy^3*invdrey);
        A(n1+2,n1+2-2)=-1/dy^4*invrey+2/2/dy^3*invdrey+(3/5)*(-1/dy^4*invrey-2/2/dy^3*invdrey);
        A(n1+2,n1+2-3)=(-1/10)*(-1/dy^4*invrey-2/2/dy^3*invdrey);
        
        %B(n1,n1+1)=i/dy^2;
        B(n1+2,n1+2+0)=-2*i/dy^2-i*k^2;
        B(n1+2,n1+2-1)=i/dy^2;
      elseif n1==Ny-1
        %A(n1,n1+2)=-1/dy^4*invrey-2/2/dy^3*invdrey;
        A(n1+2,n1+2+1)=4/dy^4*invrey+4/2/dy^3*invdrey+(i*k*Ub(n1,1)+2*k^2*invrey+1*invddrey+CU)/dy^2+(4*k^2*invdrey+dCU)/2/dy;
        A(n1+2,n1+2+0)=-6/dy^4*invrey-(i*k*Ub(n1,1)+2*k^2*invrey+1*invddrey+CU)*2/dy^2-i*k^3*Ub(n1,1)-i*k*Ub(n1,3)-k^4*invrey+k^2*invddrey;
        A(n1+2,n1+2-1)=4/dy^4*invrey-4/2/dy^3*invdrey+(i*k*Ub(n1,1)+2*k^2*invrey+1*invddrey+CU)/dy^2-(4*k^2*invdrey+dCU)/2/dy;
        A(n1+2,n1+2-2)=-1/dy^4*invrey+2/2/dy^3*invdrey;
        
        B(n1+2,n1+2+1)=i/dy^2;
        B(n1+2,n1+2+0)=-2*i/dy^2-i*k^2;
        B(n1+2,n1+2-1)=i/dy^2;
    else
        
        
        A(n1+2,1)=0;
        A(n1+2,2)=i*k*(dCU*chi+CU*dchi);
        A(n1+2,n1+2+2)=-1/dy^4*invrey-2/2/dy^3*invdrey;
        A(n1+2,n1+2+1)=4/dy^4*invrey+4/2/dy^3*invdrey+(i*k*Ub(n1,1)+2*k^2*invrey+1*invddrey+CU)/dy^2+(4*k^2*invdrey+dCU)/2/dy;
        A(n1+2,n1+2+0)=-6/dy^4*invrey-(i*k*Ub(n1,1)+2*k^2*invrey+1*invddrey+CU)*2/dy^2-i*k^3*Ub(n1,1)-i*k*Ub(n1,3)-k^4*invrey+k^2*invddrey;
        A(n1+2,n1+2-1)=4/dy^4*invrey-4/2/dy^3*invdrey+(i*k*Ub(n1,1)+2*k^2*invrey+1*invddrey+CU)/dy^2-(4*k^2*invdrey+dCU)/2/dy;
        A(n1+2,n1+2-2)=-1/dy^4*invrey+2/2/dy^3*invdrey;
        
        B(n1+2,n1+2+1)=i/dy^2;
        B(n1+2,n1+2+0)=-2*i/dy^2-i*k^2;
        B(n1+2,n1+2-1)=i/dy^2;
    end
end
[eivecs,omegavsk]=eig(full(B\A));

%[garbage,nom]=max(imag(omegavsk*ones(Ny+2,1)),[],1);
%vecmax=eivecs(:,nom);
%answer=cat(1,omegavsk(nom,nom),vecmax);

complexfreqs=cat(2,real(omegavsk*ones(Ny+2,1)),imag(omegavsk*ones(Ny+2,1)));
[garbage,indexeig]=sortrows(complexfreqs,[-2,1]);

%[garbage,indexeig2]=sortrows(eivecs',[-2]);




clear ratio
for nmode=1:Ny+2%[indexeig(1) indexeig(2) indexeig(3) indexeig(4) indexeig(5) indexeig(6)]%indexeig2(1) indexeig2(2) indexeig2(3) indexeig2(4) indexeig2(5) indexeig2(6)]
    omegar=complexfreqs(nmode,1);
    omegai=complexfreqs(nmode,2);
    qqdotdot=i*eivecs(1,nmode);
    qqdot=eivecs(2,nmode);
    qq=qqdot/(omegar+i*omegai)/i;
    epsi=eivecs(3:Ny+2,nmode)';
    for n1=1:Ny
        ev(n1,1)=-i/k*epsi(n1);
    end
    for n1=1:Ny
        if n1==1
            eu(n1,1)=i*k*(ev(n1+1,1)-0)/2/dy;
        elseif n1==Ny
            eu(n1,1)=i*k*(0-ev(n1-1,1))/2/dy;
        else
            eu(n1,1)=i*k*(ev(n1+1,1)-ev(n1-1,1))/2/dy;
        end
    end
    clear Efluid Esolid Etot
    Efluid=0;
    Esolid=0;
    for n1=1:Ny
        Efluid=Efluid+real((eu(n1,1)*conj(eu(n1,1))+ev(n1,1)*conj(ev(n1,1)))*dy);
    end
    Esolid=real(M*qqdot*conj(qqdot)+M/UR^2*qq*conj(qq));
    Etot=Esolid+Efluid;
    ratio(nmode)=Esolid/Etot;
end

[garbage,indexeig2]=sort(ratio,'descend');


nw=1;
wholeegg(nw,:)=cat(2,complexfreqs(indexeig(nw),1),complexfreqs(indexeig(nw),2),eivecs(:,indexeig(nw))',ratio(indexeig(nw)));
for nw=2:1:3;
    wholeegg(nw,:)=cat(2,complexfreqs(indexeig2(nw-1),1),complexfreqs(indexeig2(nw-1),2),eivecs(:,indexeig2(nw-1))',ratio(indexeig2(nw-1)));
end

%for nw=1:6
%    wholeegg(nw,:)=cat(2,complexfreqs(indexeig(nw),1),complexfreqs(indexeig(nw),2),eivecs(:,indexeig(nw))',ratio(indexeig(nw)));
%end
answer=wholeegg;

%wholeegg=cat(2,real(omegavsk*ones(Ny+2,1)),imag(omegavsk*ones(Ny+2,1)),eivecs');
%wholeegg=sortrows(wholeegg,[2 1]);
%answer=wholeegg(1:4,:);

%omegak=eigs(B\A,1,'li');

