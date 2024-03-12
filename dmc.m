clear all;
setenv('TZ', 'America/New_York');
fclose('all');
timeonly = false;

%what I will run
ns = [2:6];                  %uncomment parpool
ln = length(ns);
N0s = ones(1,ln) * 500;
t0s = ones(1,ln) * 50000*1;
dts = ones(1,ln) * 100*1;
C8s = ones(1,ln) * 80157;
bs = ones(1,ln) * -1;
cs = ones(1,ln) * -1;
% C8 = 80000(80157) is realistic 3sigma
% C8 = 104179       is realistic 1/4(1sigma) + 3/4(3sigma)
% C8 = 5.0000*10^12 has 0 bound states
% C8 = 2.0000*10^12 has 1 bound state
% C4 = 82.0563;
C12s = ones(1,ln) * 3.1925*10^8;
%C12 = -1 is the repulsive wall
%C12 = 3.1925*10^8  is realistic
% C6 = 1.39339*10^3
fs = ones(1,ln) * 0; %Hz
jobkeys = ["n","N0","t0","dt","C8","C12","b","c","f","datai"];

%what I will use
cores = 10;
coresPerN = 3;
runs = ln*coresPerN;
%parpool('Processes',cores);

%how I will store my results
keys = ["n","E/n","E","runtime","Emean","Rmax","N0","dE","boxsize",...
        "C8","C12","b","c","t0","dt","meanruntime","N1","f"]; %"datai"
types = length(keys);
ind = containers.Map(keys,1:types);
data = zeros(ln,types);
dataPF = cell(cores,1);

%parameters
nb = 1000;
rmin = 1.0;
rmax = 100.0; %20000;
maxextinctions = 10;
pkeys = ["cores","coresPerN","nb","rmin","rmax"];
ptypes = length(pkeys);
pind = containers.Map(pkeys,1:ptypes);
params = [cores,coresPerN,nb,rmin,rmax];

%https://www.desmos.com/calculator/b8srjezsdo
A = 0.394649;
B = 2.5886;
C = -4.16087;
D = 0.0007;

%assign job to cores
fprintf("\n\nRunning n = "+sprintf("%d, ",ns)+sprintf("\nwith %d walkers" + ...
        ", %d timesteps, C8 = %d, C12 = %d, b = %.2f, c = %.2f, f = %d.\n", ...
        N0s(1), 2*t0s(1)/dts(1), C8s(1), C12s(1), bs(1), cs(1), fs(1)));
fprintf("Assigning %d jobs to %d cores, with %d cores per n.\n\n", ...
        ln*coresPerN, cores, coresPerN);
jobs = cell(runs,1);
runtimes = zeros(runs,1);
for i = 1:ln
    n = ns(i);
    for j = 1:coresPerN
        runtimes((i-1)*coresPerN + j) = (A*n^2+B*n+C)*N0s(i) * ... 
                                        t0s(i)/dts(i)*1.0e-06 + D*n^4;
        jobs{(i-1)*coresPerN + j} = [ns(i),N0s(i),t0s(i),dts(i), ...
                                     C8s(i),C12s(i),bs(i),cs(i),fs(i),i];
        fprintf("n=%d (%.1f sec), ", ns(i), runtimes((i-1)*coresPerN+j));
        if mod((i-1)*coresPerN + j, floor(sqrt(ln*coresPerN))) == 0
            fprintf("\n");
        end
    end
end
fprintf("\n\n");
assignments = cell(cores,1);
totaltimes = zeros(cores,1);
[garbage, sortedJobIDs] = sort(runtimes, 'descend');
for i = sortedJobIDs'
    [time,core] = min(totaltimes);
    assignments{core} = [jobs{i}; assignments{core}];
    totaltimes(core) = totaltimes(core) + runtimes(i);
end
for i = 1:cores
    fprintf("Core %d will do the following jobs: ", i);
    for job = assignments{i}'
        fprintf("n=%d, ", job(1));
    end
    fprintf(". Total runtime: %.1f sec.\n", totaltimes(i));
end
[mintime, mincore] = min(totaltimes);
[maxtime, maxcore] = max(totaltimes);
fprintf("\nCore %d will run in %.1f seconds (%.2f hours), while core" + ...
        " %d will require %.1f seconds (%.2f hours).\n\n", mincore, ...
        mintime, mintime/3600, maxcore, maxtime, maxtime/3600);

%run calculations
pyrun(["import numpy", "print('Python works.\n')"]);
if timeonly
    return;
end
main = string(datetime("now","Format","M-d-y@HH.mm.ss"))+"dmc";
mkdir(main);
G = fopen(main+"/intermediate.txt", 'w');
fprintf(G,"%%[n,E/n,E,tstop,E,Rmax,N0,dE,R,C8,C12,b,c,t0,dt,tstop,N1,f]\n" +...
          "inter = [...\n");
fclose(G);

parfor corei = 1:cores
    pyrun(["import numpy",sprintf("print('On core %d, Python works.')", ...
                                  corei)]);
    H = fopen(main+"/intermediate.txt", 'a');
    ass = assignments{corei};
    [jobn, garbage] = size(ass);
    for jobi = 1:jobn
        tstart = tic;

        %initial conditions, binary search, basin hopping
        job = ass(jobi,:);
        n = job(1); N0 = job(2); t0 = job(3); dt = job(4); C8 = job(5);
        C12 = job(6); b = job(7); c = job(8); f = job(9); datai = job(10);
        Nmax = 4*N0; steps = 2*floor(t0/dt);
        sub = sprintf("n=%dcore=%djob=%d", n, corei, jobi);
        mkdir(main+"/"+sub);
        fprintf("Core=%d n=%d: Starting binary search for boxsize.\n", ...
                corei, n);
        Ra = rmin; Rb = rmax;
        N1s = zeros(maxextinctions, 1); ERs = zeros(maxextinctions, 1);
        while Rb-Ra > 0.1
            R = (Ra+Rb)/2;
            fprintf("Core=%d n=%d: Ra=%.1f, Rb=%.1f, testing R=%.1f.\n",...
                    corei, n, Ra, Rb, R);
            failed = false;
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
            ER = 0;
            extinctions = 0;
            for t = 1:maxextinctions
                for i = 1:Nmax
                    if flag(i) ~= 0
                        flag(i) = 1;
                        %1amu = 1822.89 electron masses
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
                N1s(t) = N1;
                if N1 == 0
                    if t >= 2 && N1s(t-1) > 0
                        failed = true;
                    end
                    extinctions = extinctions + 1;
                    for k = 1:N0
                        ii = q.poll();
                        flag(ii) = 1;
                        for j = 0:n-1
                            x(ii,3*j+1) = (2*rand-1)*R;
                            x(ii,3*j+2) = (2*rand-1)*R;
                            x(ii,3*j+3) = (2*rand-1)*R;
                            while x(ii,3*j+1)^2 + x(ii,3*j+2)^2 + ...
                                                  x(ii,3*j+3)^2 > R^2
                                x(ii,3*j+1) = (2*rand-1)*R;
                                x(ii,3*j+2) = (2*rand-1)*R;
                                x(ii,3*j+3) = (2*rand-1)*R;
                            end
                        end
                    end
                    ER = 0;
                    ERs(t) = ER;
                    continue;
                end
                V1 = Vsum/N1;
                ER = V1 + 1/dt*(1-N1/N0);
                ERs(t) = ER;
            end
            if (extinctions >= maxextinctions) || failed
                fprintf(sprintf("Core=%d n=%d: R=%.1f caused too " + ...
                        "many extinctions. Changing Ra=%.1f to " + ...
                        "Ra=%.1f. ", corei, n, R, Ra, R) + "N1s = " + ...
                        sprintf("%d ", N1s) + ".\n");
                Ra = R;
            elseif extinctions < maxextinctions
                fprintf(sprintf("Core=%d n=%d: R=%.1f okay. Changing" + ...
                        " Rb=%.1f to Rb=%.1f. ", corei, n, R, Rb, R) + ...
                        "N1s = " + sprintf("%d ", N1s) + ".\n");
                Rb = R;
            end
        end
        R = Rb;
        fprintf("Core=%d n=%d: Found optimal boxsize R=%.1f.\n",corei,n,R);

        %basin hopping
        fprintf("Core=%d n=%d: Matlab is calling Python to use " + ...
                "basin-hopping to optimize potential C8=%d C12=%.4e b=%.2f c=%.2f f=%d" + ...
                " for n=%d with initial R=%.1f...\n", corei,n,C8,C12,b,c,f,n,R);
        stuff = double(pyrunfile('basin.py','stuff',C8=C8,C12=C12,b=b,c=c,f=f,n=n,R=R));
        xsol = stuff(1:end-1); Vsol = stuff(end);
        fprintf("Core=%d n=%d: basin-hopping found V=%.2e.\n",corei,n,Vsol);
        attempts = 1;
        m = max(abs(xsol));
        while m > R+10 && attempts < 5
            fprintf("Core=%d n=%d: Basin-hopping optimization failed " +...
                    "for the %dth time, since one or more particles " + ...
                    "is too far away. Trying again...\n",corei,n,attempts);
            stuff = double(pyrunfile('basin.py', 'stuff', ...
                    C8=C8,C12=C12,b=b,c=c,f=f,n=n,R=R));
            attempts = attempts + 1;
            if max(abs(stuff(1:end-1))) < m
                xsol = stuff(1:end-1); Vsol = stuff(end);
                fprintf("Core=%d n=%d: basin-hopping found V=%.2e.\n",corei,n,Vsol);
                m = max(abs(xsol));
                fprintf("Core=%d n=%d: Basin-hopping optimization " + ...
                        "attempt %d was better.\n", corei, n, attempts);
            else
                fprintf("Core=%d n=%d: Basin-hopping optimization " + ...
                        "attempt %d was worse.\n", corei, n, attempts);
            end
        end
        if attempts == 5
            fprintf("Core=%d n=%d: Basin-hopping optimization failed" + ...
                    " for the 5th time. Too many attempts.\n", corei,n);
        else
            fprintf("Core=%d n=%d: Basin-hopping optimization" + ...
                    " complete.\n", corei,n);
        end
        F = fopen(main+"/"+sub+"/xyzbasin.xyz", 'w');
        fprintf(F, "%d\n\nYb %.5f %.5f %.5f\n", n,xsol(1),xsol(2),xsol(3));
        for i = 2:n
            fprintf(F, "Li %.5f %.5f %.5f\n", xsol((i-1)*3+1), ...
                    xsol((i-1)*3+2), xsol((i-1)*3+3));
        end
        fclose(F);

        %variables to update and keep track of
        ER = 0;
        ERlist = zeros(steps,1);
        ERlist(1) = ER;
        extinctions = 0;
        N1mean = 0;
        X = zeros(n,1); Y = zeros(n,1); Z = zeros(n,1);
        Vnb = 1000;
        Vcounts = zeros(Vnb,1);
        Vmin = Vsol - 0.01;
        Vmax = Vsol + 0.1;
        Rmax = 0;
        counts = zeros(nb,1);
        counted = 0; uncounted = 0;
        C4 = 82.0563;
        rmaxnow = R+5; %R+5 for high density, 20000 for low density
        r = linspace(rmin, rmaxnow, nb+1);
        r = r(2:end);
        intermin = 0;
        intermax = rmaxnow*2;
        interr = linspace(intermin, intermax, nb+1);
        interr = interr(2:end);
        interpdf = zeros(nb,1);
        intern = 10000000;
        interbins = zeros(intern,3);
        interi = 0;
        
        %rum dmc calculation
        flag = zeros(Nmax, 1);
        x = zeros(Nmax, 3*n);
        for i = 1:N0
            flag(i) = 1;
            x(i,:) = xsol;
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
            if N1 == 0
                extinctions = extinctions + 1;
                fprintf("Core=%d n=%d: for the %dth time, they all" + ...
                        " died, at timestep %d! Generating %d new" + ...
                        " walkers.\n", corei, n, extinctions, t, N0);
                if mod(extinctions,100) == 0
                    fprintf("Core=%d n=%d: they died too many times. "+ ...
                            "Increasing R=%.1f to R=%.1f.\n", corei, n, ...
                            R, R + 0.2);
                    R = R + 0.2;
                end
                for k = 1:N0
                    top = q.poll();
                    flag(top) = 1;
                    for j = 0:n-1
                        x(top,3*j+1) = (2*rand-1)*R;
                        x(top,3*j+2) = (2*rand-1)*R;
                        x(top,3*j+3) = (2*rand-1)*R;
                        while x(top,3*j+1)^2 + x(top,3*j+2)^2 + ...
                                               x(top,3*j+3)^2 > R^2
                            x(top,3*j+1) = (2*rand-1)*R;
                            x(top,3*j+2) = (2*rand-1)*R;
                            x(top,3*j+3) = (2*rand-1)*R;
                        end
                    end
                end
                ER = 0;
                continue;
            end
            N1mean = N1mean + N1;
            V1 = Vsum/N1;
            ER = V1 + 1/dt*(1-N1/N0);
            ERlist(t) = ER;
            if t >= 0.75*steps
                for i = 1:Nmax
                    if flag(i) ~= 0
                        if i == 1
                            for j = 0:n-1
                                X(j+1) = x(i,3*j+1);
                                Y(j+1) = x(i,3*j+2);
                                Z(j+1) = x(i,3*j+3);
                            end
                        end
                        Vwalk = V(x(i,:),C8,C12,b,c,f);
                        Vindex = floor((Vwalk-Vmin)/(Vmax-Vmin) * Vnb) + 1;
                        if 1 <= Vindex && Vindex <= Vnb
                            Vcounts(Vindex) = Vcounts(Vindex) + 1;
                        end
                        bins = zeros(n,1);
                        for j = 1:n-1
                            dist = sqrt((x(i,3*j+1)-x(i,1))^2 + ...
                                        (x(i,3*j+2)-x(i,2))^2 + ...
                                        (x(i,3*j+3)-x(i,3))^2);
                            bins(j) = floor((dist-rmin)/(rmaxnow-rmin)*nb) + 1;
                            Rmax = max(Rmax, dist);
                            bin = floor((dist-rmin)/(rmaxnow-rmin)*nb) + 1;
                            if 1 <= bin && bin <= nb
                                counts(bin) = counts(bin) + 1;
                                counted = counted + 1;
                            else
                                uncounted = uncounted + 1;
                            end
                        end
                        for ja = 1:n-2
                            for jb = ja+1:n-1
                                interdist = sqrt((x(i,3*ja+1)-x(i,3*jb+1))^2 ...
                                               + (x(i,3*ja+2)-x(i,3*jb+2))^2 ...
                                               + (x(i,3*ja+3)-x(i,3*jb+3))^2);
                                interbin = floor((interdist-intermin)/(intermax-intermin)*nb) + 1;
                                if interi < intern && interbin <= nb && ...
                                   bins(ja) <= nb && bins(jb) <= nb
                                    interi = interi + 1;
                                    interbins(interi,:) = [bins(ja), bins(jb), interbin];
                                end
                            end
                        end
                    end
                end
            end
            if t <= 1 || mod(t,floor(steps/10)) == 0
                fprintf("Core=%d n=%d: t = %d/%d, N1=%d, V1=%.4e, " + ...
                        "ER=%.4e.\n", corei, n, t, steps, N1, V1, ER);
            end
        end
        N1mean = floor(N1mean/steps); 
        tstop = toc(tstart);

        %store results
        E = mean(ERlist(floor(0.75*steps):end));
        dE = std(ERlist(floor(0.75*steps):end));
        fprintf("Core=%d n=%d: bins do not count %d/%d walkers, time" + ...
                " %.1f seconds, E=%.4e Hartree, E/n=%.4e Hartree/" + ...
                "particle.\n", corei, n, uncounted, uncounted+counted, ...
                tstop, E, E/n);
        % keys = ["n","E/n","E","runtime","Emean","Rmax","N0","dE", ...
        %        "boxsize","C8","C12","b","c","t0","dt","meanruntime","N1","f"]; %"datai"
        dataPF{corei} = [dataPF{corei}; [n, E/n, E, tstop, E, Rmax, N0, ...
                         dE, R, C8, C12, b, c, t0, dt, tstop, N1mean, f,  datai]];
        F = fopen(main+"/"+sub+"/data.txt",'w');
        fprintf(F, "%%DMC Yb+ ion in ultracold gas of Li atoms.\n");
        fprintf(F, "n(particles) \t%d\n", n);
        fprintf(F, "E/n(Hartree) \t%.4e\n", E/n);
        fprintf(F, "E(Hartree) \t%.4e\n", E);
        fprintf(F, "runtime(s) \t%.2f\n", tstop);
        fprintf(F, "Rmax(Bohr) \t%.8f\n", Rmax);
        fprintf(F, "N0(walkers) \t%d\n", N0);
        fprintf(F, "dE(Hartree) \t%.4e\n", dE);
        fprintf(F, "boxsize(Bohr) \t%.1f\n", R);
        fprintf(F, "C8(HBohr^8) \t%d\n", C8);
        fprintf(F, "C12(HBohr^12) \t%.4e\n", C12);
        fprintf(F, "b(Bohr) \t%.2f\n", b);
        fprintf(F, "c(Bohr) \t%.2f\n", c);
        fprintf(F, "t0(atomic time) %d\n", t0);
        fprintf(F, "dt(atomic time) %d\n", dt);
        fprintf(F, "meanruntime(s)  %.2f\n", tstop);
        fprintf(F, "N1(walkers) \t%d\n", N1mean);
        fprintf(F, "f(Hz) \t\t%d\n", f);
        %pkeys = ["cores","coresPerN","nb","rmin","rmax"];
        fprintf(F, "cores \t\t%d\n", cores);
        fprintf(F, "coresPerN \t%d\n", coresPerN);
        fprintf(F, "nb \t\t%d\n", nb);
        fprintf(F, "rmin(Bohr) \t%.1f\n", rmin);
        fprintf(F, "rmax(Bohr) \t%.1f\n", rmaxnow);
        fclose(F);
        fprintf(H, "[%d,%.4e,%.4e,%.2f,%.4e,%.8f,%d,%.4e,%.1f," + ...
                   "%d,%d,%.2f,%.2f,%d,%d,%.2f,%d,%d];\n", n, E/n, E, tstop, E, Rmax,...
                   N0, dE, R, C8, C12, b, c, t0, dt, tstop, N1mean, f);
        fseek(H, 0, 'eof');

        %plot results
        set(groot,'defaultAxesTickLabelInterpreter','latex');  
        set(groot,'defaulttextinterpreter','latex');
        set(groot,'defaultLegendInterpreter','latex');
        figure('visible','off'); clf; hold on;
        plot(ERlist);
        yline(E);
        title(["$E_R$ vs. timestep",sprintf("$n=%d$, $N_0=%d$ walkers, " + ...
              "$C_8=%d$, $C_{12}=%d$, \n$f=%d$, core$=%d$, job$=%d$, " + ...
              "$dt=%d$, \n$b=%.2f$, $c=%.2f$", ...
              n, N0, C8, C12, f, corei, jobi, dt, b, c)], ...
              "Interpreter", "latex", "FontSize", 14);
        xlabel("Timestep", "Interpreter", "latex", "FontSize", 14);
        ylabel("Reference Energy, $E_R$ (Hartree)", "Interpreter", ...
               "latex", "FontSize", 14);
        saveas(gcf, main+"/"+sub+sprintf("/energyn=%d.png",n));
        figure('visible','off'); clf; hold on;
        ax = gca;
        ax.FontSize = 25;
        ax.LineWidth = 1.2;
        box on;
        p = counts.^2' ./ r.^2;
        p = p / sum(p) / ((rmaxnow-rmin)/nb);
        plot(r, p, "LineWidth", 1.7);
        if C8 ~= -1
            xline((2*C8/C4)^(1/4), "k--", "LineWidth", 1.7);
        end
        xlabel("Atom-ion distance (Bohr)", ...
               "Interpreter", "latex", "FontSize", 25);
        ylabel(["Probability density", "$|\psi|^24\pi r^2$ (Bohr$^{-1}$)"], ...
               "Interpreter","latex","FontSize",25);
        saveas(gcf, main+"/"+sub+sprintf("/distsn=%d.png",n));
        figure('visible','off'); clf; hold on;
        ax = gca;
        ax.FontSize = 25;
        ax.LineWidth = 1.2;
        box on;
        total = 0; %for normalization later
        for i = 1:interi-1
            v = interbins(i,:);
            bina = v(1);
            binb = v(2);
            interbin = v(3);
            interpdf(interbin) = interpdf(interbin) + p(bina)*p(binb);
            total = total + p(bina)*p(binb);
        end
        dx = (intermax-intermin)/nb;
        interpdf = interpdf / total / dx;
        plot(interr, interpdf, "LineWidth", 1.7);
        C6 = 1.39339*10^3;
        if C12 ~= -1
            xline((2*C12/C6)^(1/6), "k--", "LineWidth", 1.7);
        else
            xline(30, "k--", "LineWidth", 1.7);
        end
        xlabel("Atom-atom distance (Bohr)", ...
               "Interpreter", "latex", "FontSize", 25);
        ylabel("Probability density (Bohr$^{-1}$)", ...
               "Interpreter", "latex", "FontSize", 25);
        saveas(gcf, main+"/"+sub+sprintf("/intern=%d.png",n));
        figure('visible','off');
        scatter3(X,Y,Z,100,'ko');
        hold on; scatter3(X(1),Y(1),Z(1),400,'r.'); hold off;
        title(["Positions of first walker's Yb$^+$ ion and Li atoms", ...
              "during last timestep at equilibrium", ...
              sprintf("$n=%d$, $N_0=%d$ walkers, $C_8=%d$, $C_{12}=%d$, " + ...
              "\n$f=%d$, core$=%d$, job$=%d$, $dt=%d$, \n$b=%.2f$, $c=%.2f$", ...
              n, N0, C8, C12, f, corei, jobi, dt, b, c)], "Interpreter", ...
              "latex", "FontSize", 14);
        saveas(gcf, main+"/"+sub+sprintf("/geometryn=%d.png",n));
        F = fopen(main+"/"+sub+sprintf("/xyzn=%dcore=%djob=%d.xyz", ...
                  n, corei, jobi), 'w');
        fprintf(F, "%d\n\nYb %.6f %.6f %.6f\n", n, X(1), Y(1), Z(1));
        for j = 2:n
            fprintf(F, "Li %.6f %.6f %.6f\n", X(j), Y(j), Z(j));
        end
        fclose(F);
    end
    fprintf("Core %d finished!\n",corei);
    fclose(H);
end
G = fopen(main+"/intermediate.txt",'a');
fprintf(G, "];");
fclose(G);

%organize results
%keys = ["n","E/n","E","runtime","Emean","Rmax","N0","dE","boxsize", ...
%        "C8","C12","b","c","t0","dt","meanruntime","N1","f"]; %"datai"
for i = 1:ln
    data(i,:)=[ns(i),100,100,0,0,0,N0s(i),0,rmax,C8s(i),C12s(i),bs(i),cs(i),...
               t0s(i),dts(i),0,0,fs(i)];
end
energies = cell(ln,1);
for corei = 1:cores
    for d = dataPF{corei}'
        i = d(end);
        data(i,ind("E/n")) = min(data(i,ind("E/n")),d(ind("E/n")));
        update = false;
        if d(ind("E")) < data(i,ind("E"))
            data(i,ind("E")) = min(data(i,ind("E")),d(ind("E")));
            update = true;
        end
        data(i,ind("runtime")) = max(data(i,ind("runtime")), ...
                                     d(ind("runtime")));
        data(i,ind("Emean")) = data(i,ind("Emean")) + d(ind("Emean"));
        data(i,ind("Rmax")) = max(data(i,ind("Rmax")), d(ind("Rmax")));
        if update
            data(i,ind("dE")) = d(ind("dE"));
        end
        data(i,ind("boxsize")) = min(data(i,ind("boxsize")), ...
                                 d(ind("boxsize")));
        energies{i} = [energies{i}, d(ind("E"))];
        data(i,ind("meanruntime")) = data(i,ind("meanruntime")) + ...
                                     d(ind("meanruntime"));
        data(i,ind("N1")) = data(i,ind("N1")) + d(ind("N1"));
    end
end
for i = 1:ln
    data(i,ind("Emean")) = data(i,ind("Emean")) / coresPerN;
    data(i,ind("meanruntime")) = data(i,ind("meanruntime")) / coresPerN;
    data(i,ind("N1")) = floor(data(i,ind("N1")) / coresPerN);
end

%write results to data files
save(main+"/work.mat");
F = fopen(main+"/data.txt", 'w');
fprintf(F, "%%DMC Yb+ ion in ultracold gas of Li atoms\n");
%keys = ["n","E/n","E","runtime","Emean","Rmax","N0","dE","boxsize", ...
%        "C8","C12","b","c","t0","dt","meanruntime","N1","f"]; %"datai"
fprintf(F, "%%[n,E/n,E,runtime,Emean,Rmax,N0,dE,boxsize,C8,C12,b,c,t0,dt," + ...
           "meanruntime,N1,f]\n");
s = "data{100} = [...\n";
for i = 1:ln
    s = s + "[";
    for j = 1:types
        s = s + string(data(i,j)) + ",";
    end
    s = s + "];\n";
end
s = s + "];\n";
fprintf(F, s);
%pkeys = ["cores","coresPerN","nb","rmin","rmax"];
fprintf(F, "%%[cores,coresPerN,nb,rmin,rmax]\n");
s = "params{100} = [";
for i = 1:ptypes
    s = s + string(params(i)) + ",";
end
s = s + "];\n";
fprintf(F, s);
fclose(F);
F = fopen(main+"/energies.txt", 'a');
for i = 1:ln
    fprintf(F, "n=%d, N0=%d, t=%d, dt=%d, C8=%d, C12=%.4e, b=%.2f, c=%.2f, f=%d\n", ...
            ns(i), N0s(i), t0s(i), dts(i), C8s(i), C12s(i), bs(i), cs(i), fs(i));
    e = energies{i};
    for j = 1:coresPerN
        fprintf(F, "%.4e, ", e(j));
    end
    fprintf(F, "\n");
end
fclose(F);

%play sound when finished
[garbage1, garbage2] = audioread('sound.mp3');
player = audioplayer(garbage1, garbage2);
play(player);

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