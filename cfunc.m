function answer=cfunc(y,H)


%DRAG TERM PROFILE. dataC IS EQUIVALENT TO THE PROFILE OF FIG 2(b)
%Profile from Sylvain Dupont 18 mars 2009
datay=[0,5.797101557E-002,1.739130467E-001,2.898550630E-001,4.057970941E-001,5.217391253E-001,6.376811862E-001,7.536231875E-001,8.695651889E-001,9.855072498E-001,1.101449251E+000,1.217391253E+000];
dataC=[4.75,4.750000000E+000,4.750000000E+000,4.750000000E+000,4.750000000E+000,4.750000000E+000,4.750000000E+000,4.300000191E+000,3.349999905E+000,2.400000095E+000,1.450000048E+000,0.000000000E+000];

CD=0.2*2;
h=0.69;
dataC=dataC*CD*h;

hdata=max(datay);

if y<hdata
    answer(1,1)=pchip(datay,dataC,y);
    delta=(1.739130467E-001-5.797101557E-002)/h;
    uc=answer(1,1);
    up=pchip(datay,dataC,y+delta);
    if y<delta
        um=pchip(datay,dataC,delta);
    else
        um=pchip(datay,dataC,y-delta);
    end
    answer(1,2)=(up-um)/2/delta;
else
    answer(1,1)=0;
    answer(1,2)=0;
end

 return
%}
%TANH Profile
%{
constant=0.2*3;
%constant=1*.01*0.69/.05^2;
delta=0.4;

answer(1,1)=constant*(.5-.5*tanh(2*(y-1)/delta));
answer(1,2)=constant*(-1.0*(1-tanh(2*(y-1)/delta)^2)/delta);
%}