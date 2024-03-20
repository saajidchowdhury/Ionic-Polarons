clear all;
main = string(datetime("now","Format","M-d-y@HH.mm.ss"))+"dmcSimple";
mkdir(main);

n = 3;
N0 = 500;
t0 = 50000;
dt = 100;

Nmax = 4*N0; 
steps = 2*t0/dt;
R = 10;

ER = 0;
ERlist = zeros(steps,1);
ERlist(1) = ER;

%rum dmc calculation
flag = zeros(Nmax, 1);
x = zeros(Nmax, 3*n);
for i = 1:N0
    flag(i) = 1;
    for j = 0:n-1
        x(i,3*j+1) = (2*rand-1)*R;
        x(i,3*j+2) = (2*rand-1)*R;
        x(i,3*j+3) = (2*rand-1)*R;
        while x(i,3*j+1)^2 + x(i,3*j+2)^2 + x(i,3*j+3)^2 > R^2
            x(i,3*j+1) = (2*rand-1)*R;
            x(i,3*j+2) = (2*rand-1)*R;
            x(i,3*j+3) = (2*rand-1)*R;
        end
    end
end
q = java.util.LinkedList();
for i = N0+1:Nmax
    q.add(i);
end
for t = 1:steps
    for i = 1:Nmax
        if flag(i) ~= 0
            flag(i) = 1;
            mLi = 7.0160034366 * 1822.89;
            mYb = 173.938859 * 1822.89;
            for j = 1:3
                x(i,j) = x(i,j) + sqrt(dt/mYb)*randn;
            end
            for j = 4:3*n
                x(i,j) = x(i,j) + sqrt(dt/mLi)*randn;
            end
        end
    end
    N1 = 0;
    Vsum = 0;
    for i = 1:Nmax
        if flag(i) == 1
            W = exp(-(V(x(i,:))-ER)*dt);
            mn = min(floor(W+rand),3);
            if mn == 0
                flag(i) = 0;
                q.add(i);
            elseif mn == 2
                top = q.poll();
                flag(top) = 2;
                x(top,:) = x(i,:);
            elseif mn == 3
                top1 = q.poll();
                top2 = q.poll();
                flag(top1) = 2;
                x(top1,:) = x(i,:);
                flag(top2) = 2;
                x(top2,:) = x(i,:);
            end
            N1 = N1 + mn;
            Vsum = Vsum + mn*V(x(i,:));
        end
    end
    V1 = Vsum/N1;
    ER = V1 + 1/dt*(1-N1/N0);
    ERlist(t) = ER;
    if t <= 1 || mod(t,floor(steps/10)) == 0
        fprintf("n=%d: t = %d/%d, N1=%d, V1=%.4e, " + ...
                "ER=%.4e.\n", n, t, steps, N1, V1, ER);
    end
end

E = mean(ERlist(floor(0.75*steps):end));
dE = std(ERlist(floor(0.75*steps):end));
fprintf("n=%d: E=%.4e Hartree, E/n=%.4e Hartree/particle.\n", n, E, E/n);

%plot results
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
figure('visible','on'); clf; hold on;
plot(ERlist, "LineWidth", 2);
yline(E, "LineWidth", 2);
xlabel("Timestep", "FontSize", 20);
ylabel("Reference Energy, $E_R$ (Hartree)", "FontSize", 20);
ax = gca; ax.FontSize = 20;
saveas(gcf, main+"/energy.png");

%potential
function V = V(x)
    n = length(x)/3;
    V = 0;
    C4 = 82.0563; C8 = 80157;
    for j = 1:n-1
        R = sqrt((x(3*j+1)-x(1))^2+(x(3*j+2)-x(2))^2+(x(3*j+3)-x(3))^2);
        V = V - C4/R^4 + C8/R^8;
    end
    C6 = 1.39339*10^3; C12 = 3.1925*10^8;
    for ja = 1:n-2
        for jb = ja+1:n-1
            R = sqrt((x(3*ja+1)-x(3*jb+1))^2 + (x(3*ja+2)-x(3*jb+2))^2 ...
                   + (x(3*ja+3)-x(3*jb+3))^2);
            V = V - C6/R^6 + C12/R^12;
        end
    end
end
