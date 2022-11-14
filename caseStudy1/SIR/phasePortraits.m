[I,S] = generatePhaseDomain(1,1,50,1);

sPrime = zeros(size(I));
iPrime = zeros(size(I));

betaParam = 0.1;
rho = 0.4;

differentialVector = @(t,Y,betaParam,rho) ...
    [-Y(1)*Y(2)*betaParam; betaParam*Y(2)*(Y(1)-rho)];

t=0; 
for i = 1:numel(S)
    Yprime = differentialVector(t,[S(i); I(i)],betaParam,rho);
    sPrime(i) = Yprime(1);
    iPrime(i) = Yprime(2);
end

quiver(S,I,sPrime,iPrime,'Color','r','AutoScaleFactor',2);
hold on;
xline(rho,'b--');
%figure(gcf)
xlabel('S')
ylabel('I')
formatSpec = 'Phase portrait w/ rho=%.1f, beta=%.1f';
title(sprintf(formatSpec,rho,betaParam))
legend('phase portrait','S = \rho')
axis tight equal;