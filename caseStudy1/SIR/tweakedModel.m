N = 10000;
betaParam = 0.3*N;
v = 1*N;
vaccineParam = v/(betaParam*N^2);
gammaParam = 0.1;
muParam = gammaParam/(betaParam*N);
R = 1/(muParam*N)
f = @(t,x) [-x(1)*x(2)-vaccineParam;x(1)*x(2)-x(2)*muParam;x(2)*muParam+vaccineParam];
g = @(t,x) [-x(1)*x(2);x(1)*x(2)-x(2)*muParam;x(2)*muParam];

[t,xa]=ode45(f,[0 6], [99 1 0]);
[T,ya]=ode45(g,[0 100*betaParam*N], [0.999 1-0.999 0]);

figure(1)
plot(T,ya(:,1),'m')
hold on;
plot(T,ya(:,2),'c')
plot(T,ya(:,3),'k')
legend('S','I','R')
hold off
%%
figure(2)
plot(t,xa(:,1),'r')
hold on
plot(t,xa(:,2),'g')
plot(t,xa(:,3),'b')
hold off