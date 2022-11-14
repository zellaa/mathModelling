function [y1_out,y2_out,y3_out] = odsolver(t_vals,y_init,f_1,f_2,f_3)
%Author: Zella Baig, Date: 30/12/20
%This function solves an ODE using the RK4 method
%inputs:
%t_vals: this is the array of time values
%y_init: the 3x1 array of initial y values
%f_1,f_2,f_3: the anonymous functions which govern the lorentz equations
%output: this outputs 3 Nx1 vectors, corresponding to values of y_i for
%each given time in t_vals
%
%Usage: [y1_out,y2_out,y3_out] = odsolver(t,y_initial_values,function1,function2,function3)

%first work out the time interval
dt = t_vals(2)-t_vals(1);
f1=f_1;
f2=f_2;
f3=f_3;

%initialise variables
t=t_vals;  
y1(1)=y_init(1);
y2(1)=y_init(2);
y3(1)=y_init(3);

%create a loop to perform the RK4 method for each time interval
for i=1:(length(t)-1) 
    
    %find fi_0:
    f1_0=f1(t(i),y1(i),y2(i),y3(i));
    f2_0=f2(t(i),y1(i),y2(i),y3(i));
    f3_0=f3(t(i),y1(i),y2(i),y3(i));
    
    %find fi_1
    f1_1=f1(t(i)+dt,(y1(i)+0.5*f1_0*dt),(y2(i)+(0.5*f2_0*dt)),(y3(i)+(0.5*f3_0*dt)));     
    f2_1=f2(t(i)+dt,(y1(i)+0.5*f1_0*dt),(y2(i)+(0.5*f2_0*dt)),(y3(i)+(0.5*f3_0*dt)));
    f3_1=f3(t(i)+dt,(y1(i)+0.5*f1_0*dt),(y2(i)+(0.5*f2_0*dt)),(y3(i)+(0.5*f3_0*dt)));
     
    %find fi_2
    f1_2=f1(t(i)+dt,(y1(i)+0.5*f1_1*dt),(y2(i)+(0.5*f2_1*dt)),(y3(i)+(0.5*f3_1*dt)));
    f2_2=f2(t(i)+dt,(y1(i)+0.5*f1_1*dt),(y2(i)+(0.5*f2_1*dt)),(y3(i)+(0.5*f3_1*dt)));
    f3_2=f3(t(i)+dt,(y1(i)+0.5*f1_1*dt),(y2(i)+(0.5*f2_1*dt)),(y3(i)+(0.5*f3_1*dt)));
    
    %find fi_3
    f1_3=f1(t(i)+dt,(y1(i)+f1_2*dt),(y2(i)+f2_2*dt),(y3(i)+f3_2*dt));
    f2_3=f2(t(i)+dt,(y1(i)+f1_2*dt),(y2(i)+f2_2*dt),(y3(i)+f3_2*dt));
    f3_3=f3(t(i)+dt,(y1(i)+f1_2*dt),(y2(i)+f2_2*dt),(y3(i)+f3_2*dt));
    
    %find the proceding yi values
    y1(i+1) = y1(i) + dt*(f1_0 +2*f1_1 +2*f1_2 +f1_3)/6;
    y2(i+1) = y2(i) + dt*(f2_0 +2*f2_1 +2*f2_2 +f2_3)/6;
    y3(i+1) = y3(i) + dt*(f3_0 +2*f3_1 +2*f3_2 +f3_3)/6;
end

%store outputs
y1_out = y1;
y2_out = y2;
y3_out = y3;
end

