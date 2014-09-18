% Plotting eigenvectors

R = load('R.txt');
n = length(R);
rho_max = 6;
rho_min = 0;
h =(rho_max-rho_min)/(n+1);
i = linspace(0, n-1, n);
rho = rho_min + i*h;

% Read from program:
R3 = R(:,33).^2;
R7 = R(:,141).^2;
R11 = R(:,93).^2;


plot(rho, R3, 'r');
xlabel('index');
ylabel('Value');
hold on
plot(rho, R7, 'm');
plot(rho, R11);
title('Eigenvectors vs rho');
legend('Energy: 3','Energy : 7','Energy: 11');
show()

hold off
