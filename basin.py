from math import sqrt
from math import pi
from numpy.random import rand
from numpy import append
from scipy.optimize import basinhopping

#matlab arguments: C8, C12, b, c, f, n, R

def V(x, C8, C12, b, c, f):
    n = len(x)//3
    ans = 0;
    C4 = 82.0563
    for j in range(1,n,1):
        R = sqrt((x[3*j+0]-x[0])**2+(x[3*j+1]-x[1])**2+(x[3*j+2]-x[2])**2)
        if C8 == -1:
            ans = ans - C4*(R*R-c*c)/(R*R+c*c)*1/((b*b+R*R)**2)
        elif b == -1 and c == -1:
            ans = ans - C4/(R**4) + C8/(R**8)
        else:
            print("Error -1")
    C6 = 1.39339*(10**3)
    for ja in range(1,n-1,1):
        for jb in range(ja+1,n,1):
            R = sqrt((x[3*ja+0]-x[3*jb+0])**2+(x[3*ja+1]-x[3*jb+1])**2+(x[3*ja+2]-x[3*jb+2])**2)
            if C12 == -1:
                if R < 30:
                    ans = ans + 0.015
            else:
                ans = ans - C6/(R**6) + C12/(R**12)
    m = 173.938859 * 1822.89
    w = 2*pi*f * 2.4188843265857*(10**(-17))
    R = sqrt(x[0]**2+x[1]**2+x[2]**2)
    ans = ans + 0.5*m*(w**2)*(R**2)
    return ans

n = int(n)
x = [0 for i in range(3*n)]
for j in range(0,n,1):
    x[3*j+0] = (2*rand()-1)*R
    x[3*j+1] = (2*rand()-1)*R
    x[3*j+2] = (2*rand()-1)*R

res = basinhopping(lambda z : V(z, C8, C12, b, c, f), x, stepsize=0.1*R, niter=1)

xsol = res.x
val = V(xsol, C8, C12, b, c, f)
stuff = append(xsol, val)
