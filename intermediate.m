%[n,E/n,E,tstop,E,Rmax,N0,dE,R,C8,C12,b,c,t0,dt,tstop,N1,f]
inter = [...
[2,-1.0226e-02,-2.0452e-02,6.29,-2.0452e-02,8.19914826,500,9.9657e-05,2.7,80157,319250000,-1.00,-1.00,50000,100,6.29,499,0];
[2,-1.0224e-02,-2.0447e-02,6.83,-2.0447e-02,8.05586401,500,1.0106e-04,2.8,80157,319250000,-1.00,-1.00,50000,100,6.83,499,0];
[2,-1.0224e-02,-2.0448e-02,6.06,-2.0448e-02,8.06549321,500,1.0041e-04,2.7,80157,319250000,-1.00,-1.00,50000,100,6.06,499,0];
[3,-1.4086e-02,-4.2258e-02,7.96,-4.2258e-02,8.17588438,500,1.2901e-04,3.8,80157,319250000,-1.00,-1.00,50000,100,7.96,497,0];
[5,-1.7746e-02,-8.8730e-02,10.73,-8.8730e-02,8.03529474,500,1.6921e-04,5.9,80157,319250000,-1.00,-1.00,50000,100,10.73,493,0];
[5,-1.7746e-02,-8.8728e-02,11.11,-8.8728e-02,8.03809666,500,1.7848e-04,5.7,80157,319250000,-1.00,-1.00,50000,100,11.11,493,0];
[5,-1.7753e-02,-8.8764e-02,11.00,-8.8764e-02,8.22707691,500,1.9212e-04,5.9,80157,319250000,-1.00,-1.00,50000,100,11.00,493,0];
[6,-1.8869e-02,-1.1321e-01,13.50,-1.1321e-01,8.06812521,500,1.9975e-04,6.3,80157,319250000,-1.00,-1.00,50000,100,13.50,491,0];
[6,-1.8861e-02,-1.1317e-01,13.94,-1.1317e-01,8.09614295,500,2.1272e-04,6.8,80157,319250000,-1.00,-1.00,50000,100,13.94,490,0];
[6,-1.8863e-02,-1.1318e-01,13.81,-1.1318e-01,8.06808791,500,1.9989e-04,6.4,80157,319250000,-1.00,-1.00,50000,100,13.81,490,0];
[3,-1.4072e-02,-4.2216e-02,5.99,-4.2216e-02,8.05526120,500,1.3322e-04,3.9,80157,319250000,-1.00,-1.00,50000,100,5.99,497,0];
[4,-1.6340e-02,-6.5361e-02,6.26,-6.5361e-02,8.12001134,500,1.6906e-04,4.8,80157,319250000,-1.00,-1.00,50000,100,6.26,495,0];
[4,-1.6341e-02,-6.5364e-02,5.99,-6.5364e-02,8.11779501,500,1.5583e-04,4.6,80157,319250000,-1.00,-1.00,50000,100,5.99,495,0];
[4,-1.6336e-02,-6.5343e-02,6.88,-6.5343e-02,8.07913902,500,1.6198e-04,4.7,80157,319250000,-1.00,-1.00,50000,100,6.88,495,0];
[3,-1.4079e-02,-4.2236e-02,3.49,-4.2236e-02,8.03833399,500,1.1525e-04,3.8,80157,319250000,-1.00,-1.00,50000,100,3.49,497,0];
];

[s, i] = sort(inter(:,1));
inter = inter(i,:);
disp(inter);
ns = unique(s);
ln = length(ns);
ls = length(s);

keys = ["n","E/n","E","runtime","Emean","Rmax","N0","dE","boxsize", ...
        "C8","C12","b","c","t0","dt","meanruntime","N1","f"]; %"datai"
types = length(keys);
ind = containers.Map(keys, 1:types);

energies = cell(ln,1);
data = zeros(ln,types);
for i = 1:ln
    data(i,:) = [ns(i),100,100,0,0,0,0,0,100,0,0,0,0,0,0,0,0,0];
end
i = 0;
for j = 1:ls
    if j == 1 || s(j) ~= s(j-1)
        i = i+1;
    end
    energies{i} = [energies{i}, inter(j,ind("E"))];
    data(i,ind("n")) = inter(j,ind("n"));
    data(i,ind("E/n")) = min(data(i,ind("E/n")), inter(j,ind("E/n")));
    data(i,ind("E")) = min(data(i,ind("E")), inter(j,ind("E")));
    data(i,ind("runtime")) = max(data(i,ind("runtime")), inter(j,ind("runtime")));
    update = false;
    if inter(j,ind("Emean")) < data(i,ind("Emean"))
        update = true;
        data(i,ind("Emean")) = inter(j,ind("Emean"));
    end
    data(i,ind("Rmax")) = max(data(i,ind("Rmax")), inter(j,ind("Rmax")));
    data(i,ind("N0")) = inter(j,ind("N0"));
    if update
        data(i,ind("dE")) = inter(j,ind("dE"));
    end
    data(i,ind("boxsize")) = min(data(i,ind("boxsize")), inter(j,ind("boxsize")));
    data(i,ind("C8")) = inter(j,ind("C8"));
    data(i,ind("C12")) = inter(j,ind("C12"));
    data(i,ind("b")) = inter(j,ind("b"));
    data(i,ind("c")) = inter(j,ind("c"));
    data(i,ind("t0")) = inter(j,ind("t0"));
    data(i,ind("dt")) = inter(j,ind("dt"));
    data(i,ind("meanruntime")) = data(i,ind("meanruntime")) + inter(j,ind("meanruntime"));
    data(i,ind("N1")) = data(i,ind("N1")) + inter(j,ind("N1"));
    data(i,ind("f")) = inter(j,ind("f"));
end
for i = 1:ln
    l = length(energies{i});
    data(i,ind("Emean")) = data(i,ind("Emean")) / l;
    data(i,ind("meanruntime")) = data(i,ind("meanruntime")) / l;
    data(i,ind("N1")) = data(i,ind("N1")) / l;
end

fprintf("\n\n\n%%DMC Yb+ ion in ultracold gas of Li atoms\n");
%keys = ["n","E/n","E","runtime","Emean","Rmax","N0","dE","boxsize", ...
%        "C8","C12","b","c","t0","dt","meanruntime","N1","f"]; %"datai"
fprintf("%%[n,E/n,E,runtime,Emean,Rmax,N0,dE,boxsize,C8,C12,b,c,t0,dt,meanruntime,N1,f]\n");
s = "data{100} = [...\n";
for i = 1:ln
    s = s + "[";
    for j = 1:types
        s = s + string(data(i,j)) + ",";
    end
    s = s + "];\n";
end
s = s + "];\n";
fprintf(s);
