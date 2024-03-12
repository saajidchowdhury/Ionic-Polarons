clear all;
main = string(datetime("now","Format","M-d-y@HH.mm.ss"))+"dmcSimple";
mkdir(main);

n = 3;
N0 = 500;
t0 = 50000;
dt = 100;
C4 = 82.0563;
C8 = 80157;
b = -1;
c = -1;
C6 = 1.39339*10^3;
C12 = 3.1925*10^8;
f = 0; %Hz

Nmax = 4*N0; 
steps = 2*floor(t0/dt);
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
            %1amu = 1822.89 electron mass
            %mass of Li7 = 7.0160034366(45) amu
            %mass of Yb174 = 173.938859 amu
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
            W = exp(-(V(x(i,:),C8,C12,b,c,f)-ER)*dt);
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
            Vsum = Vsum + mn*V(x(i,:),C8,C12,b,c,f);
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




function V = V(x, C8, C12, b, c, f)
    n = length(x)/3;
    V = 0;
    C4 = 82.0563;
    for j = 1:n-1
        R = sqrt((x(3*j+1)-x(1))^2+(x(3*j+2)-x(2))^2+(x(3*j+3)-x(3))^2);
        if C8 == -1
            V = V - C4*(R^2-c^2)/(R^2+c^2)*1/(b^2+R^2)^2;
        elseif b == -1 && c == -1
            V = V - C4/R^4 + C8/R^8;
        else
            fprintf("Error -1.\n");
        end
    end
    C6 = 1.39339*10^3;
    for ja = 1:n-2
        for jb = ja+1:n-1
            R = sqrt((x(3*ja+1)-x(3*jb+1))^2 + (x(3*ja+2)-x(3*jb+2))^2 ...
                   + (x(3*ja+3)-x(3*jb+3))^2);
            if C12 == -1
                if R < 30
                    V = V + 0.015;
                end
            else
                V = V - C6/R^6 + C12/R^12;
            end
        end
    end

    %trap
    m = 173.938859 * 1822.89;
    w = 2*pi*f * 2.4188843265857*10^(-17);
    R = sqrt(x(1)^2+x(2)^2+x(3)^2);
    V = V + 0.5*m*w^2*R^2;
end