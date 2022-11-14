%Author: Zella Baig, Date: 30/12/20
%This script generates the final y data across time, and also allows
%plotting of various things.
%No inputs required, but ensure dt, numsteps, a, b, r, and y_init are
%defined
%output: this an Nx3 array containing the yi values against time, and also
%plots should the user wish.
%
%Usage: simply run the script after (un)commenting whichever plots needed.

%find the final t value, i.e. the max value to check in the ODE
numsteps = 1000;
%the time interval
dt=0.02;
t=transpose(0:dt:dt*numsteps);

%define the value of the constants
a=10;   
b=8/3;
r=28;

%create the initial value(s)
y_init = [4;5;6];
y_init2 = [4.01,5.01,6.01];

%first column is y1, second is y2, third is y3
y = solve_lorentz(y_init,a,b,r,t);
y2 = solve_lorentz(y_init2,a,b,r,t);

%create plots

%uncomment this section to view y3 against y2
figure(1)
plot(y(:,2),y(:,3));

%display info about the parameters used
%avals = sprintf('a = %.f',a);
%bvals = sprintf('b = %.3f',b);
%rvals = sprintf('r = %.f',r);
%info_vals = {avals,bvals,rvals};
%text(5,5,info_vals)

%title('y_3 against y_2');
%xlabel('y_2 values');
%ylabel('y_3 values');


%figure(2)
%used to display comparison plot, comment out if needed
%subplot(2,1,1)
%plot values
%plot(t,y(:,1),'m');
%hold on;
%plot(t,y2(:,1),'b')

title('y_1 against t, r=28 comparison');
xlabel('t values');
ylabel('y_1 values');
legend('y_0 = (4,5,6)','y_0 = (4.01,5.01,6.01)','Location', 'NorthEast');

%now do the same for the y values with a measurement error
%subplot(2,1,2)
%plot(t,y2(:,1),'b');

%uncomment for parameter information
%avals = sprintf('a = %.f',a);
%bvals = sprintf('b = %.3f',b);
%rvals = sprintf('r = %.f',r);
%info_vals = {avals,bvals,rvals};
%text(5,5,info_vals)

%title('y_1 against t, r=28 comparison');
%xlabel('t values');
%ylabel('y_1 values');
%legend('y_0 = (4.01,5.01,6.01)','y_0 = (4.01,5.01,6.01)','Location', 'NorthEast');
