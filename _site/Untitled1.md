
## Hartree Fock calculation of H2O system in Python
1. This tutorial shows how hartree fock calculation is done typically. This does not intend to be the most efficient algorithmic implementation of HF. But a pedagogical tutorial on what is Hartree Fock and what is the mathematical equations behind it.
2. We also visualize the Atomic orbitla and Molecular Orbitals, showing how the linear combination of atomic orbitals (LCAO) gave molecular orbitals (MO)
3. We implemented the McMurchie-Davidson algorithm (using recursion) to calculate the 1-electron and 2-electron integrals (this integrals are saved in Memory, I didn't try to save it on the disk and read them later)
4. Scipy linear algebra package is used to solve the Hartree Fock Roothaan-Hall equation (generalized eigenvalue-eigenvector problem) (FC=SCe)

### Import libraries


```python
import numpy as np
import collections
import copy
from scipy.special import comb
from scipy.special import factorial2
from scipy.special import factorial
```

### Define Atomic Orbital object, Molecular Orbital object, and Gaussian Primitive object. Give two functions that convert between "S, P, D, F" to "(0,0,0),(1,0,0), (2,0,0)...."



```python
def getAngularQuanta(orbitaltype):
    if orbitaltype == "S":
        return 0, 0, 0
    if orbitaltype == "Px":
        return 1, 0, 0
    if orbitaltype == "Py":
        return 0, 1, 0
    if orbitaltype == "Pz":
        return 0, 0, 1
    if orbitaltype == "Dx2":
        return 2, 0, 0
    if orbitaltype == "Dy2":
        return 0, 2, 0
    if orbitaltype == "Dz2":
        return 0, 0, 2
    if orbitaltype == "Dxy":
        return 1, 1, 0
    if orbitaltype == "Dyz":
        return 0, 1, 1
    if orbitaltype == "Dzx":
        return 1, 0, 1
    
def getorbitalType(angular):
    if angular == (0,0,0):
        return "S"
    if angular == (1,0,0):
        return "Px"
    if angular == (0,1,0):
        return "Py"
    if angular == (0,0,1):
        return "Pz"
    if angular == (2,0,0):
        return "Dx2"
    if angular == (0,2,0):
        return "Dy2"
    if angular == (0,0,2):
        return "Dz2"
    if angular == (1,1,0):
        return "Dxy"
    if angular == (0,1,1):
        return "Dyz"
    if angular == (1,0,1):
        return "Dzx"
    
class Gprimitive: #gaussian primitive class for only one variable. The total will be product of gprimitive(x)*gprimitive(y)*gprimitive(z)
    def __init__(self, angular, center, exponent):
        self.angular = angular
        self.center = center
        self.exponent = exponent
    
    def __call__(self, x):
        return (x - self.center)**self.angular*np.exp(-exponent*(x-self.center)**2)
    
    def __repr__(self):
        return str(self.center)+str(self.angular)+str(self.exponent)
        
    
class Ao(object): #atomic orbital
    def __init__(self, center, angular, contract_num):
        self.center = center # the center of the atomic orbital
        self.exponents = [0 for i in range(contract_num)] #list of gaussian primitive exponents
        self.coeffs = [0 for i in range(contract_num)] #list of gaussian primivite coeffs
        self.angular = angular #angular momentum could be S, Px, Py, Pz, Dx2, Dy2, Dz2, Dxy, Dyz, Dzx, ... (0,0,0), (1,0,0),(0,1,0)...
        
    def __repr__(self):
        return getorbitalType(self.angular) + str(self.center) + str(self.exponents) + str(self.coeffs)
    
    def __call__(self, x, y, z):
        res = 0
        x0, y0, z0 = self.center
        l, m, n = self.angular
        for i in range(len(self.coeffs)):
            exponent = self.exponents[i]
            gprimitivex = Gprimitive(l, x0, exponent)
            gprimitivey = Gprimitive(m, y0, exponent)
            gprimitivez = Gprimitive(n, z0, exponent)
            res += self.coeffs[i]*gprimitivex(x)*gprimitivey(y)*gprimitivez(z)
        return res
            
    
        
class Mo(object): #molecular orbital->linear combination of atomic orbitals
    def __init__(self, aolist, coeffs):
        self.aolist = aolist
        self.coeffs = coeffs
        
    def __repr__(self):
        return str(self.coeffs) + repr(self.aolist)
    
    def __call__(self, x, y, z):
        res = 0
        for i in range(len(self.aolist)):
            ao = aolist[i]
            res += self.coeffs[i]*ao(x, y, z)
        return res
    
    
```

# Read Basis Set Information from Basis set data. The data is in Gaussian 94 Format and downloaded from EMSL basis set exchange
Gaussian Format and Basis Set downloaded from EMSL Basis Set Exchange. In this tutorial 6-311++G(2d,2p) is used.


```python
file = open("basis.dat")
aodict = collections.defaultdict(list)
newatomtype = False
contract_num = None
origin = (0,0,0)
orbitalType = None
for line in file:
    if line.strip() == "****":
        newatomtype = True
        linecount = 0
    elif line.strip() != "****" and newatomtype == True:
        if linecount == 0:
            atomtype, _ = line.split()
            linecount += 1
            print(atomtype)
        elif linecount == 1:
            orbitalType, contract_num, _ = line.split()
            contract_num = int(contract_num)
            if orbitalType == "S":
                aos = Ao(origin, (0,0,0), contract_num)
                aodict[atomtype].append(aos)
            elif orbitalType == "P":
                aopx = Ao(origin,(1,0,0), contract_num)
                aopy = Ao(origin,(0,1,0), contract_num)
                aopz = Ao(origin,(0,0,1), contract_num)
                aodict[atomtype].append(aopx)
                aodict[atomtype].append(aopy)
                aodict[atomtype].append(aopz)
            elif orbitalType == "D":
                aodx2 = Ao(origin,(2,0,0), contract_num)
                aody2 = Ao(origin,(0,2,0), contract_num)
                aodz2 = Ao(origin,(0,0,2), contract_num)
                aodxy = Ao(origin,(1,1,0), contract_num)
                aodyz = Ao(origin,(0,1,1), contract_num)
                aodzx = Ao(origin,(1,0,1), contract_num)
                aodict[atomtype].append(aodx2)
                aodict[atomtype].append(aody2)
                aodict[atomtype].append(aodz2)
                aodict[atomtype].append(aodxy)
                aodict[atomtype].append(aodyz)
                aodict[atomtype].append(aodzx)
            elif orbitalType == "SP":
                aos = Ao(origin,(0,0,0), contract_num)
                aodict[atomtype].append(aos)
                aopx = Ao(origin,(1,0,0), contract_num)
                aopy = Ao(origin,(0,1,0), contract_num)
                aopz = Ao(origin,(0,0,1), contract_num)
                aodict[atomtype].append(aopx)
                aodict[atomtype].append(aopy)
                aodict[atomtype].append(aopz)
            linecount += 1
            print(orbitalType, contract_num)
        elif contract_num and 1 <linecount <= 1 + contract_num:
            if orbitalType == "S" or orbitalType == "P" or orbitalType == "D":
                exponent, coeff = line.split()
                exponent = float(exponent)
                coeff = float(coeff)
                if orbitalType == "S":
                    aos.exponents[linecount - 2] = exponent
                    aos.coeffs[linecount - 2] = coeff
                elif orbitalType == "P":
                    aopx.exponents[linecount - 2] = exponent
                    aopy.exponents[linecount - 2] = exponent
                    aopz.exponents[linecount - 2] = exponent
                    aopx.coeffs[linecount - 2] = coeff
                    aopy.coeffs[linecount - 2] = coeff
                    aopz.coeffs[linecount - 2] = coeff
                elif orbitalType == "D":
                    aodx2.exponents[linecount - 2] = exponent
                    aody2.exponents[linecount - 2] = exponent
                    aodz2.exponents[linecount - 2] = exponent
                    aodxy.exponents[linecount - 2] = exponent
                    aodyz.exponents[linecount - 2] = exponent
                    aodzx.exponents[linecount - 2] = exponent
                    aodx2.coeffs[linecount - 2] = coeff
                    aody2.coeffs[linecount - 2] = coeff
                    aodz2.coeffs[linecount - 2] = coeff
                    aodxy.coeffs[linecount - 2] = coeff
                    aodyz.coeffs[linecount - 2] = coeff
                    aodzx.coeffs[linecount - 2] = coeff
                print(exponent, coeff)
            elif orbitalType == "SP":
                exponent, coeffs, coeffp = line.split()
                exponent = float(exponent)
                coeffs = float(coeffs)
                coeffp = float(coeffp)
                aos.exponents[linecount - 2] = exponent
                aos.coeffs[linecount - 2] = coeffs
                aopx.exponents[linecount - 2] = exponent
                aopy.exponents[linecount - 2] = exponent
                aopz.exponents[linecount - 2] = exponent
                aopx.coeffs[linecount - 2] = coeffp
                aopy.coeffs[linecount - 2] = coeffp
                aopz.coeffs[linecount - 2] = coeffp
                print(exponent, coeffs, coeffp)
            if linecount < 1 + int(contract_num):
                linecount += 1
            else:
                linecount = 1
                
print(aodict)
```

    H
    S 2
    5.447178 0.156285
    0.824547 0.904691
    S 1
    0.183192 1.0
    O
    S 3
    322.037 0.0592394
    48.4308 0.3515
    10.4206 0.707658
    SP 2
    7.40294 -0.404453 0.244586
    1.5762 1.22156 0.853955
    SP 1
    0.373684 1.0 1.0
    defaultdict(<class 'list'>, {'H': [S(0, 0, 0)[5.447178, 0.824547][0.156285, 0.904691], S(0, 0, 0)[0.183192][1.0]], 'O': [S(0, 0, 0)[322.037, 48.4308, 10.4206][0.0592394, 0.3515, 0.707658], S(0, 0, 0)[7.40294, 1.5762][-0.404453, 1.22156], Px(0, 0, 0)[7.40294, 1.5762][0.244586, 0.853955], Py(0, 0, 0)[7.40294, 1.5762][0.244586, 0.853955], Pz(0, 0, 0)[7.40294, 1.5762][0.244586, 0.853955], S(0, 0, 0)[0.373684][1.0], Px(0, 0, 0)[0.373684][1.0], Py(0, 0, 0)[0.373684][1.0], Pz(0, 0, 0)[0.373684][1.0]]})


## Read Geometry information
Read Geometry information in Cartesian Coordinate. The unit is in Angstrom and the full molecule is used. No symmetry is used. In this tutorial H2O is used as illustration.


```python
geomfile = open("geom.dat")
aolist = []
atomlist = []
for line in geomfile:
    atomtype, x, y, z = line.split()
    x = float(x)
    y = float(y)
    z = float(z)
    x = x * 1.8897261339213 # from Angstrom to atomic unit
    y = y * 1.8897261339213 # from Angstrom to atomic unit
    z = z * 1.8897261339213 # ....
    atomlist.append((atomtype, (x, y, z)))
    if atomtype in aodict:
        for ao in aodict[atomtype]:
            aonew = copy.copy(ao)
            aonew.center = [x,y,z]
            aolist.append(aonew)
            
print(aolist)
print(len(aolist))
```

    [S[-0.0963760328299863, -0.5801459231138391, -0.4875493425516954][322.037, 48.4308, 10.4206][0.0592394, 0.3515, 0.707658], S[-0.0963760328299863, -0.5801459231138391, -0.4875493425516954][7.40294, 1.5762][-0.404453, 1.22156], Px[-0.0963760328299863, -0.5801459231138391, -0.4875493425516954][7.40294, 1.5762][0.244586, 0.853955], Py[-0.0963760328299863, -0.5801459231138391, -0.4875493425516954][7.40294, 1.5762][0.244586, 0.853955], Pz[-0.0963760328299863, -0.5801459231138391, -0.4875493425516954][7.40294, 1.5762][0.244586, 0.853955], S[-0.0963760328299863, -0.5801459231138391, -0.4875493425516954][0.373684][1.0], Px[-0.0963760328299863, -0.5801459231138391, -0.4875493425516954][0.373684][1.0], Py[-0.0963760328299863, -0.5801459231138391, -0.4875493425516954][0.373684][1.0], Pz[-0.0963760328299863, -0.5801459231138391, -0.4875493425516954][0.373684][1.0], S[-0.7540007274345988, 1.1565123939598356, -0.9297452578892796][5.447178, 0.824547][0.156285, 0.904691], S[-0.7540007274345988, 1.1565123939598356, -0.9297452578892796][0.183192][1.0], S[-0.6557349684706911, -0.5763664708459965, -1.417294600440975][5.447178, 0.824547][0.156285, 0.904691], S[-0.6557349684706911, -0.5763664708459965, -1.417294600440975][0.183192][1.0]]
    13


### Define a Function that do implicit plot of 3D surfaces to visualize the Atomic and Molecular orbtials.


```python
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np

def plot_implicit(fn, isovalue, bbox=(-2.5,2.5),elev=0, azim=30):
    ''' create a plot of an implicit function
    fn  ...implicit function (plot where fn==0)
    bbox ..the x,y,and z limits of plotted interval'''
    xmin, xmax, ymin, ymax, zmin, zmax = bbox*3
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    A = np.linspace(xmin, xmax, 100) # resolution of the contour
    B = np.linspace(xmin, xmax, 15) # number of slices
    A1,A2 = np.meshgrid(A,A) # grid on which the contour is plotted

    for z in B: # plot contours in the XY plane
        X,Y = A1,A2
        Z = fn(X,Y,z) - isovalue
        cset = ax.contour(X, Y, Z+z, [z], zdir='z', colors='red')
        # [z] defines the only level to plot for this contour for this value of z

    for y in B: # plot contours in the XZ plane
        X,Z = A1,A2
        Y = fn(X,y,Z) - isovalue
        cset = ax.contour(X, Y+y, Z, [y], zdir='y', colors='red')

    for x in B: # plot contours in the YZ plane
        Y,Z = A1,A2
        X = fn(x,Y,Z) - isovalue
        cset = ax.contour(X+x, Y, Z, [x], zdir='x', colors='red')
        
## now plot the negative part
    for z in B: # plot contours in the XY plane
        X,Y = A1,A2
        Z = fn(X,Y,z) + isovalue
        cset = ax.contour(X, Y, Z+z, [z], zdir='z', colors='blue')
        # [z] defines the only level to plot for this contour for this value of z

    for y in B: # plot contours in the XZ plane
        X,Z = A1,A2
        Y = fn(X,y,Z) + isovalue
        cset = ax.contour(X, Y+y, Z, [y], zdir='y', colors='blue')

    for x in B: # plot contours in the YZ plane
        Y,Z = A1,A2
        X = fn(x,Y,Z) + isovalue
        cset = ax.contour(X+x, Y, Z, [x], zdir='x', colors='blue')

    # must set plot limits because the contour will likely extend
    # way beyond the displayed level.  Otherwise matplotlib extends the plot limits
    # to encompass all values in the contour.
    ax.set_zlim3d(zmin,zmax)
    ax.set_xlim3d(xmin,xmax)
    ax.set_ylim3d(ymin,ymax)
    ax.view_init(elev, azim)
    plt.show()
```

## Now Start visualize some of the atomic orbitals


```python
for i in range(len(aolist)):
    ao = aolist[i]
    print(i, getorbitalType(ao.angular))
plot_implicit(aolist[4], 0.05,[-4,4],0,90) #Pz orbital on Oxygen.
```

    0 S
    1 S
    2 Px
    3 Py
    4 Pz
    5 S
    6 Px
    7 Py
    8 Pz
    9 S
    10 S
    11 S
    12 S



![png](output_12_1.png)


### The major work horse of the SCF.  In this section, we define the function that calculates 1-electron and 2-electron integrals using Boys Function and McMurchie-Davidson recurrence formulation.


```python
def gint(m, center1, exponent1, n, center2, exponent2): #calculate one electron integral of two gaussian primitives
    newcenter = (exponent1 * center1 + exponent2 * center2) / (exponent1 + exponent2)
    newexponent = exponent1 + exponent2
    tempexponent = exponent1 * exponent2 / (exponent1 + exponent2)
    e12 = np.exp(-tempexponent*(center1 - center2)**2)
    res = 0
    for i in range(m + 1):
        for j in range(n + 1):
            if (i + j) % 2 == 0:
                res += (np.sqrt(np.pi/newexponent)*comb(m, i)*comb(n, j)*
                        factorial2(i + j - 1)/(2*newexponent)**((i + j)/2)*
                        (newcenter - center1)**(m - i)*(newcenter - center2)**(n - j)
                       )
    res = e12 * res
    return res
                        

def E(l1, l2, t, center1, center2, exponent1, exponent2):#calculate the gaussian-hermite expansion coefficient using recurence
    newcenter = (exponent1 * center1 + exponent2 * center2) / (exponent1 + exponent2)
    sumexponent = exponent1 + exponent2
    diffcenter = center1 - center2
    redexponent = exponent1 * exponent2 / (exponent1 + exponent2)
    if t > l1 + l2:
        return 0
    if l1 < 0 or l2 < 0 or t < 0:
        return 0
    elif l1 == 0 and l2 == 0 and t == 0:
        return np.exp(-redexponent*diffcenter**2)
    elif l1 > 0:
        return (1/(2*sumexponent)*E(l1-1, l2, t - 1,center1,center2,exponent1,exponent2)
                +(newcenter - center1)*E(l1-1,l2,t,center1,center2,exponent1,exponent2)
                +(t + 1)*E(l1-1,l2, t+1, center1,center2,exponent1,exponent2))
    elif l1 == 0:
        return (1/(2*sumexponent)*E(l1, l2-1, t - 1,center1,center2,exponent1,exponent2)
                +(newcenter - center2)*E(l1,l2-1,t,center1,center2,exponent1,exponent2)
                +(t + 1)*E(l1,l2-1, t+1, center1,center2,exponent1,exponent2))
    return 0


def S(m1, m2, center1, center2, exponent1, exponent2): #calculate overlap type integral
    return np.sqrt(np.pi/(exponent1 + exponent2))*E(m1, m2, 0, center1, center2, exponent1, exponent2)


def T(m1, m2, center1, center2, exponent1, exponent2): #calculate kinetic type integral
    res = 0
    res += -2*exponent2*S(m1, m2 + 2, center1, center2, exponent1, exponent2)
    res += exponent2*(2*m2+1)*S(m1, m2, center1, center2, exponent1, exponent2)
    res += -1/2*m2*(m2-1)*S(m1, m2 - 2, center1, center2, exponent1, exponent2)
    return res
    

def F(n, x): #calculate Boys function value by numerical integration
    if x < 1e-7:
        return 1/(2*n + 1)
    if n == 20:
        res1 = 1/(2*n + 1)
        #if x < 1e-7:
        #    return res1
        for k in range(1,11):
            res1 += (-x)**k/factorial(k)/(2*n+2*k+1)
        res2 = factorial2(2*n-1)/2**(n+1)*np.sqrt(np.pi/x**(2*n+1))
        res = min(res1, res2)
        return res
    return (2*x*F(n+1,x)+np.exp(-x))/(2*n+1)


def R(t, u, v, n, p, x, y, z):
    if t < 0 or u < 0 or v < 0:
        return 0
    if t == 0 and u == 0 and v == 0:
        return (-2*p)**n*F(n,p*(x**2+y**2+z**2))
    if t > 0:
        return (t-1)*R(t-2,u,v,n+1,p,x,y,z)+x*R(t-1,u,v,n+1,p,x,y,z)
    if u > 0:
        return (u-1)*R(t,u-2,v,n+1,p,x,y,z)+y*R(t,u-1,v,n+1,p,x,y,z)
    if v > 0:
        return (v-1)*R(t,u,v-2,n+1,p,x,y,z)+z*R(t,u,v-1,n+1,p,x,y,z)
    
    
def overlap(ao1, ao2): #calculate overlap matrix <psi|phi>
    l1, m1, n1 = ao1.angular
    l2, m2, n2 = ao2.angular
    x1, y1, z1 = ao1.center
    x2, y2, z2 = ao2.center
    res = 0
    for i in range(len(ao1.coeffs)):
        for j in range(len(ao2.coeffs)):
            exponent1 = ao1.exponents[i]
            exponent2 = ao2.exponents[j]
            #res += (ao1.coeffs[i]*ao2.coeffs[j]*
            #        gint(l1, x1, exponent1, l2, x2, exponent2)*
            #        gint(m1, y1, exponent1, m2, y2, exponent2)*
            #        gint(n1, z1, exponent1, n2, z2, exponent2))
            res += (ao1.coeffs[i]*ao2.coeffs[j]*
                    S(l1, l2, x1, x2, exponent1, exponent2)*
                    S(m1, m2, y1, y2, exponent1, exponent2)*
                    S(n1, n2, z1, z2, exponent1, exponent2))
    return res




def kinetic(ao1, ao2): #calculate kinetic integral <psi|-1/2*del^2|phi>
    l1, m1, n1 = ao1.angular
    l2, m2, n2 = ao2.angular
    x1, y1, z1 = ao1.center
    x2, y2, z2 = ao2.center
    res = 0
    for i in range(len(ao1.coeffs)):
        for j in range(len(ao2.coeffs)):
            exponent1 = ao1.exponents[i]
            exponent2 = ao2.exponents[j]
            res += (ao1.coeffs[i]*ao2.coeffs[j]*
                    (T(l1,l2,x1,x2,exponent1,exponent2)*S(m1,m2,y1,y2,exponent1,exponent2)*S(n1,n2,z1,z2,exponent1,exponent2) +
                     S(l1,l2,x1,x2,exponent1,exponent2)*T(m1,m2,y1,y2,exponent1,exponent2)*S(n1,n2,z1,z2,exponent1,exponent2) +
                     S(l1,l2,x1,x2,exponent1,exponent2)*S(m1,m2,y1,y2,exponent1,exponent2)*T(n1,n2,z1,z2,exponent1,exponent2))
                   )
    return res


def oneelectron(ao1, centerC, ao2):
    l1, m1, n1 = ao1.angular
    l2, m2, n2 = ao2.angular
    a = l1 + m1 + n1
    b = l2 + m2 + n2
    c = a + b
    x1, y1, z1 = ao1.center
    x2, y2, z2 = ao2.center
    xc, yc, zc = centerC # coordinate of atom with charge Z
    res = 0
    for i in range(len(ao1.coeffs)):
        for j in range(len(ao2.coeffs)):
            exponent1 = ao1.exponents[i]
            exponent2 = ao2.exponents[j]
            p = exponent1 + exponent2
            xp = (exponent1*x1+exponent2*x2)/p
            yp = (exponent1*y1+exponent2*y2)/p
            zp = (exponent1*z1+exponent2*z2)/p
            xpc = xp - xc
            ypc = yp - yc
            zpc = zp - zc
            for t in range(c+1):
                for u in range(c+1):
                    for v in range(c+1):
                        res += (ao1.coeffs[i]*ao2.coeffs[j]*
                                2*np.pi/p*E(l1,l2,t,x1,x2,exponent1,exponent2)*
                                E(m1,m2,u,y1,y2,exponent1,exponent2)*
                                E(n1,n2,v,z1,z2,exponent1,exponent2)*
                                R(t,u,v,0,p,xpc,ypc,zpc))
    return res
                  
    
def twoelectron(ao1, ao2, ao3, ao4):
    res = 0
    l1, m1, n1 = ao1.angular
    l2, m2, n2 = ao2.angular
    l3, m3, n3 = ao3.angular
    l4, m4, n4 = ao4.angular
    x1, y1, z1 = ao1.center
    x2, y2, z2 = ao2.center
    x3, y3, z3 = ao3.center
    x4, y4, z4 = ao4.center
    a = l1 + m1 + n1
    b = l2 + m2 + n2
    c = l3 + m3 + n3
    d = l4 + m4 + n4
    for i in range(len(ao1.coeffs)):
        for j in range(len(ao2.coeffs)):
            for k in range(len(ao3.coeffs)):
                for l in range(len(ao4.coeffs)):
                    exponent1 = ao1.exponents[i]
                    exponent2 = ao2.exponents[j]
                    exponent3 = ao3.exponents[k]
                    exponent4 = ao4.exponents[l]
                    p = (exponent1 + exponent2)
                    q = (exponent3 + exponent4)
                    alpha = p*q/(p + q)
                    xp = (x1*exponent1+x2*exponent2)/p
                    yp = (y1*exponent1+y2*exponent2)/p
                    zp = (z1*exponent1+z2*exponent2)/p
                    xq = (x3*exponent3+x4*exponent4)/q
                    yq = (y3*exponent3+y4*exponent4)/q
                    zq = (z3*exponent3+z4*exponent4)/q
                    xpq = xp - xq
                    ypq = yp - yq
                    zpq = zp - zq
                    for t in range(a + b + 1):
                        for u in range(a + b + 1):
                            for v in range(a + b + 1):
                                for tau in range(c + d + 1):
                                    for miu in range(c + d + 1):
                                        for phi in range(c + d + 1):
                                            res += (ao1.coeffs[i]*ao2.coeffs[j]*ao3.coeffs[k]*ao4.coeffs[l]*
                                                    2*np.pi**(5/2)/p/q/np.sqrt(p+q)*
                                                    E(l1, l2, t, x1, x2, exponent1, exponent2)*
                                                    E(m1, m2, u, y1, y2, exponent1, exponent2)*
                                                    E(n1, n2, v, z1, z2, exponent1, exponent2)*
                                                    E(l3, l4, tau, x3, x4, exponent3, exponent4)*
                                                    E(m3, m4, miu, y3, y4, exponent3, exponent4)*
                                                    E(n3, n4, phi, z3, z4, exponent3, exponent4)*
                                                    (-1)**(tau + miu + phi)*
                                                    R(t+tau, u+miu, v+phi, 0, alpha, xpq, ypq, zpq)
                                                    )
    return res
                            
    
    
def atomtype2Charge(atomtype):
    if atomtype == "O":
        return 8
    if atomtype == "H":
        return 1
    if atomtype == "C":
        return 6
    if atomtype == "N":
        return 7
    if atomtype == "S":
        return 16
    
def coulombicAttraction(ao1, atomlist, ao2):
    res = 0
    for i in range(len(atomlist)):
        atomtype, centerC = atomlist[i]
        Z = atomtype2Charge(atomtype)
        res += (-Z)*oneelectron(ao1, centerC, ao2)
    return res


```

## Now Build S Matrix using the overlap integral


```python
n = len(aolist)
Smatrix = np.zeros((n,n))
for i in range(n):
    for j in range(n):
        Smatrix[i][j] = overlap(aolist[i], aolist[j])
    
print(Smatrix)
```

    [[ 3.62762722e-02  9.95722893e-02  1.93673973e-19  5.37206593e-21
       4.47877187e-19  1.16908729e-01  1.54276710e-18  1.23358108e-17
       5.84923568e-18  6.10913341e-03  6.22165339e-02  4.05905413e-02
       9.70264391e-02]
     [ 9.95722893e-02  1.29603405e+00  1.42684721e-17 -1.07089068e-18
       0.00000000e+00  2.39434572e+00  0.00000000e+00  0.00000000e+00
      -1.38677717e-16  2.24970293e-01  1.54657288e+00  8.51936046e-01
       2.31527171e+00]
     [ 1.93673973e-19  1.42684721e-17  1.20081274e-01  1.33198625e-34
       0.00000000e+00  0.00000000e+00  4.51862346e-01  0.00000000e+00
       0.00000000e+00 -3.67968817e-02 -7.72600583e-02 -1.22819876e-01
      -9.85352674e-02]
     [ 5.37206593e-21 -1.07089068e-18  1.33198625e-34  1.20081274e-01
       0.00000000e+00  0.00000000e+00  0.00000000e+00  4.51862346e-01
       0.00000000e+00  9.71733744e-02  2.04028717e-01  8.29864027e-04
       6.65778834e-04]
     [ 4.47877187e-19  0.00000000e+00  0.00000000e+00  0.00000000e+00
       1.20081274e-01 -9.69453238e-17  0.00000000e+00  0.00000000e+00
       4.51862346e-01 -2.47427308e-02 -5.19507288e-02 -2.04146551e-01
      -1.63781593e-01]
     [ 1.16908729e-01  2.39434572e+00  0.00000000e+00  0.00000000e+00
      -9.69453238e-17  8.61832884e+00 -1.19603339e-16  0.00000000e+00
       4.78413356e-16  1.52207840e+00  8.56140029e+00  2.87854472e+00
       1.15939955e+01]
     [ 1.54276710e-18  0.00000000e+00  4.51862346e-01  0.00000000e+00
       0.00000000e+00 -1.19603339e-16  5.76578663e+00  0.00000000e+00
      -6.63931903e-33 -6.91617246e-01 -1.85212767e+00 -1.11368404e+00
      -2.13339716e+00]
     [ 1.23358108e-17  0.00000000e+00  0.00000000e+00  4.51862346e-01
       0.00000000e+00  0.00000000e+00  0.00000000e+00  5.76578663e+00
       0.00000000e+00  1.82642600e+00  4.89110726e+00  7.52489217e-03
       1.44148456e-02]
     [ 5.84923568e-18 -1.38677717e-16  0.00000000e+00  0.00000000e+00
       4.51862346e-01  4.78413356e-16 -6.63931903e-33  0.00000000e+00
       5.76578663e+00 -4.65052976e-01 -1.24539619e+00 -1.85112347e+00
      -3.54605203e+00]
     [ 6.10913341e-03  2.24970293e-01 -3.67968817e-02  9.71733744e-02
      -2.47427308e-02  1.52207840e+00 -6.91617246e-01  1.82642600e+00
      -4.65052976e-01  2.25610819e+00  5.04483605e+00  5.73298304e-01
       3.09593009e+00]
     [ 6.22165339e-02  1.54657288e+00 -7.72600583e-02  2.04028717e-01
      -5.19507288e-02  8.56140029e+00 -1.85212767e+00  4.89110726e+00
      -1.24539619e+00  5.04483605e+00  2.51084590e+01  3.09593009e+00
       1.86434852e+01]
     [ 4.05905413e-02  8.51936046e-01 -1.22819876e-01  8.29864027e-04
      -2.04146551e-01  2.87854472e+00 -1.11368404e+00  7.52489217e-03
      -1.85112347e+00  5.73298304e-01  3.09593009e+00  2.25610819e+00
       5.04483605e+00]
     [ 9.70264391e-02  2.31527171e+00 -9.85352674e-02  6.65778834e-04
      -1.63781593e-01  1.15939955e+01 -2.13339716e+00  1.44148456e-02
      -3.54605203e+00  3.09593009e+00  1.86434852e+01  5.04483605e+00
       2.51084590e+01]]


## Now build T Matrix using the kinetic integral


```python
Tmatrix = np.zeros((n,n))
for i in range(n):
    for j in range(n):
        Tmatrix[i][j] = kinetic(aolist[i], aolist[j])
print(Tmatrix)
```

    [[ 1.56594600e+00  4.78962978e-02  6.78908861e-18  4.22064249e-20
       2.27868157e-17  1.19389138e-01  2.61562630e-18  2.09132274e-17
       9.87340471e-18 -1.78867914e-02 -4.92373005e-02  2.50637951e-02
       8.01241472e-03]
     [ 3.43169953e+00  2.67653865e+00  9.87769972e-17  1.70224326e-16
       0.00000000e+00  1.26286244e+00  0.00000000e+00  0.00000000e+00
      -1.26224305e-16 -2.40376917e-01 -1.29261541e+00  6.87680289e-01
      -2.67661699e-01]
     [ 7.28176557e-18  1.74603270e-17  7.22034552e-01  2.52382758e-33
       0.00000000e+00  0.00000000e+00  4.14183613e-01  0.00000000e+00
       0.00000000e+00 -8.77758511e-03 -7.40780850e-02 -2.85476643e-01
      -1.66050624e-01]
     [ 5.16328248e-18 -1.17093707e-16  2.52382758e-33  7.22034552e-01
       0.00000000e+00  0.00000000e+00  0.00000000e+00  4.14183613e-01
       0.00000000e+00  2.31798871e-02  1.95625747e-01  1.92889623e-03
       1.12196368e-03]
     [-1.74215738e-17  0.00000000e+00  0.00000000e+00  0.00000000e+00
       7.22034552e-01 -1.57856939e-17  0.00000000e+00  0.00000000e+00
       4.14183613e-01 -5.90216930e-03 -4.98111261e-02 -4.74508474e-01
      -2.76003065e-01]
     [ 4.02359929e+00  3.74486998e+00  0.00000000e+00  0.00000000e+00
      -3.72194533e-16 -3.26589848e+00  7.55390768e-17  0.00000000e+00
      -3.02156307e-16 -3.19847345e-02 -8.89098584e+00  1.17885250e+00
      -7.32217219e+00]
     [ 4.15219877e-17  0.00000000e+00  1.84953755e+00  0.00000000e+00
       0.00000000e+00  1.64926785e-16 -3.64155552e+00  0.00000000e+00
       1.08325801e-32 -4.49357986e-01 -5.62261276e-01 -1.21442179e+00
      -1.51582970e+00]
     [ 3.26095844e-16  0.00000000e+00  0.00000000e+00  1.84953755e+00
       0.00000000e+00  0.00000000e+00  0.00000000e+00 -3.64155552e+00
       0.00000000e+00  1.18666664e+00  1.48482216e+00  8.20555263e-03
       1.02420926e-02]
     [ 1.18328784e-16 -9.52467827e-17  0.00000000e+00  0.00000000e+00
       1.84953755e+00 -6.59707140e-16  1.08325801e-32  0.00000000e+00
      -3.64155552e+00 -3.02154508e-01 -3.78072238e-01 -2.01856595e+00
      -2.51955477e+00]
     [ 2.04383059e-01  2.11458960e-01 -1.25157853e-01  3.30517433e-01
      -8.41578665e-02 -1.68671733e+00  6.76615664e-01 -1.78680976e+00
       4.54965705e-01  2.94423220e+00  5.04572120e-02 -1.51881931e-01
      -2.45033733e+00]
     [ 2.12607292e+00  1.98995307e+00 -2.68742014e-01  7.09695146e-01
      -1.80705837e-01 -1.01604601e+01  3.29952851e+00 -8.71341006e+00
       2.21864848e+00  9.70000310e-01 -2.38636820e+01  6.94714612e-03
      -2.32695311e+01]
     [ 1.39374740e+00  1.35463291e+00 -5.47657338e-01  3.70038742e-03
      -9.10295305e-01 -6.49531399e-01  1.00773990e-01 -6.80905340e-04
       1.67502714e-01 -1.51881931e-01 -2.45033733e+00  2.94423220e+00
       5.04572120e-02]
     [ 3.32042071e+00  3.12403406e+00 -3.51732442e-01  2.37657055e-03
      -5.84636356e-01 -1.14464781e+01  3.37499366e+00 -2.28040112e-02
       5.60978677e+00  6.94714612e-03 -2.32695311e+01  9.70000310e-01
      -2.38636820e+01]]


## Now build 1-Electron Coulombic Attraction Matrix using 1-electron integral


```python
coulombAttractionMatrix = np.zeros((n,n))
for i in range(n):
    for j in range(n):
        coulombAttractionMatrix[i][j] = coulombicAttraction(aolist[i], atomlist, aolist[j])
print(coulombAttractionMatrix)
```

    [[-1.77274472e+00 -3.36133464e+00  3.46788189e-03 -9.67761921e-04
       5.40934821e-03 -3.93525332e+00  4.96477720e-03 -1.33033155e-03
       7.76499514e-03 -1.96486391e-01 -2.07083435e+00 -1.36749327e+00
      -3.23813099e+00]
     [-3.36133464e+00 -2.12738044e+01  9.79437480e-02 -8.20013141e-02
       1.32234912e-01 -3.29727352e+01  3.14504890e-01 -2.58026151e-01
       4.26602808e-01 -2.64177831e+00 -1.96403707e+01 -1.22258217e+01
      -2.97687359e+01]
     [ 3.46788189e-03  9.79437480e-02 -1.55033337e+00  5.00966826e-03
      -1.71420855e-02  2.28346096e-01 -4.56891809e+00  2.51850126e-02
      -7.02873816e-02  3.68067594e-01  8.70865432e-01  1.46549031e+00
       1.15552947e+00]
     [-9.67761921e-04 -8.20013141e-02  5.00966826e-03 -1.55196470e+00
       3.43247294e-03 -1.82717143e-01  2.51850126e-02 -4.58683522e+00
       1.71920073e-02 -9.33927806e-01 -1.98838526e+00 -8.60270673e-02
      -1.35433954e-01]
     [ 5.40934821e-03  1.32234912e-01 -1.71420855e-02  3.43247294e-03
      -1.56615190e+00  3.11471674e-01 -7.02873816e-02  1.71920073e-02
      -4.63155872e+00  2.61797394e-01  7.02591124e-01  2.40727839e+00
       1.87272154e+00]
     [-3.93525332e+00 -3.29727352e+01  2.28346096e-01 -1.82717143e-01
       3.11471674e-01 -7.81590971e+01  1.41774279e+00 -1.24170553e+00
       1.89354472e+00 -1.16705504e+01 -6.52346379e+01 -2.82900302e+01
      -9.12349322e+01]
     [ 4.96477720e-03  3.14504890e-01 -4.56891809e+00  2.51850126e-02
      -7.02873816e-02  1.41774279e+00 -3.60903412e+01  3.33633657e-01
      -5.59513719e-01  4.67632572e+00  1.16785306e+01  8.73784488e+00
       1.38063811e+01]
     [-1.33033155e-03 -2.58026151e-01  2.51850126e-02 -4.58683522e+00
       1.71920073e-02 -1.24170553e+00  3.33633657e-01 -3.65549043e+01
       2.26251190e-01 -1.20567202e+01 -2.83217459e+01 -4.11247894e-01
      -2.22796124e+00]
     [ 7.76499514e-03  4.26602808e-01 -7.02873816e-02  1.71920073e-02
      -4.63155872e+00  1.89354472e+00 -5.59513719e-01  2.26251190e-01
      -3.65253384e+01  3.25434920e+00  8.79931320e+00  1.43913722e+01
       2.21463363e+01]
     [-1.96486391e-01 -2.64177831e+00  3.68067594e-01 -9.33927806e-01
       2.61797394e-01 -1.16705504e+01  4.67632572e+00 -1.20567202e+01
       3.25434920e+00 -1.42921174e+01 -2.97382773e+01 -4.74408391e+00
      -2.02964730e+01]
     [-2.07083435e+00 -1.96403707e+01  8.70865432e-01 -1.98838526e+00
       7.02591124e-01 -6.52346379e+01  1.16785306e+01 -2.83217459e+01
       8.79931320e+00 -2.97382773e+01 -1.23833572e+02 -2.57019412e+01
      -1.08193177e+02]
     [-1.36749327e+00 -1.22258217e+01  1.46549031e+00 -8.60270673e-02
       2.40727839e+00 -2.82900302e+01  8.73784488e+00 -4.11247894e-01
       1.43913722e+01 -4.74408391e+00 -2.57019412e+01 -2.13191742e+01
      -4.13050760e+01]
     [-3.23813099e+00 -2.97687359e+01  1.15552947e+00 -1.35433954e-01
       1.87272154e+00 -9.12349322e+01  1.38063811e+01 -2.22796124e+00
       2.21463363e+01 -2.02964730e+01 -1.08193177e+02 -4.13050760e+01
      -1.49144815e+02]]


## Now Build 2-electron Coulombic Repulsion Tensor
### This Section takes the longest time. In the early days, the 2-electron integral is calculated and saved on disk. They will be read in when neccessary. Nowadays, with the fast computing power, it is actually faster to calculate the 2-electron integrals than read from disk.
### Here we calculate them and save them in a rank-4 tensor in memory (O(N<sup>4</sup>). The 2-electron integral make the Hartree Fock Algorithm to be O(N<sup>4</sup>) time complexity. But in reality, most of the 2-electron integrals are Zero. It can be optimized to be O(N<sup>2</sup>)


```python
import time, sys
from IPython.display import clear_output
def update_progress(progress):
    bar_length = 20
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
    if progress < 0:
        progress = 0
    if progress >= 1:
        progress = 1
    block = int(round(bar_length * progress))
    clear_output(wait = True)
    text = "Progress: [{0}] {1:.1f}%".format( "#" * block + "-" * (bar_length - block), progress * 100)
    print(text)

    
coulombRepulsionTensor = np.zeros((n,n,n,n))
count = 0
for i in range(n):
    for j in range(n):
        for k in range(n):
            for l in range(n):
                coulombRepulsionTensor[i][j][k][l] = twoelectron(aolist[i],aolist[j],aolist[k],aolist[l])
                count += 1
                update_progress(count / n**4)
                
update_progress(1)
print(coulombRepulsionTensor)
```

    Progress: [####################] 100.0%
    [[[[ 5.25865931e-03  1.14911426e-02  2.22135775e-20 ...  7.07882741e-03
         4.64671403e-03  1.10598789e-02]
       [ 1.14911426e-02  8.20034204e-02  6.88977021e-19 ...  7.67107360e-02
         4.49946748e-02  1.16458608e-01]
       [ 2.22135775e-20  6.88977021e-19  5.94422277e-03 ... -2.77337017e-03
        -4.89744402e-03 -3.57050084e-03]
       ...
       [ 7.07882741e-03  7.67107360e-02 -2.77337017e-03 ...  4.27110568e-01
         9.08275712e-02  3.83671838e-01]
       [ 4.64671403e-03  4.49946748e-02 -4.89744402e-03 ...  9.08275712e-02
         7.15552870e-02  1.46801848e-01]
       [ 1.10598789e-02  1.16458608e-01 -3.57050084e-03 ...  3.83671838e-01
         1.46801848e-01  5.40067365e-01]]
    
      [[ 1.14911426e-02  2.70372687e-02  3.16023786e-20 ...  1.66912063e-02
         1.09393747e-02  2.60657172e-02]
       [ 2.70372687e-02  2.15524108e-01  1.67989285e-18 ...  2.04720415e-01
         1.19661680e-01  3.10532067e-01]
       [ 3.16023786e-20  1.67989285e-18  1.60557730e-02 ... -7.57262516e-03
        -1.33448242e-02 -9.74746816e-03]
       ...
       [ 1.66912063e-02  2.04720415e-01 -7.57262516e-03 ...  1.16872315e+00
         2.46923793e-01  1.04751554e+00]
       [ 1.09393747e-02  1.19661680e-01 -1.33448242e-02 ...  2.46923793e-01
         1.94158375e-01  3.99191612e-01]
       [ 2.60657172e-02  3.10532067e-01 -9.74746816e-03 ...  1.04751554e+00
         3.99191612e-01  1.47362662e+00]]
    
      [[ 2.22135775e-20  3.16023786e-20  1.74187480e-04 ... -2.72777201e-05
        -6.67469497e-05 -3.62691325e-05]
       [ 3.16023786e-20  5.55697522e-19  1.84477971e-03 ... -5.52598476e-04
        -1.13376852e-03 -7.21640607e-04]
       [ 1.74187480e-04  1.84477971e-03  6.56129339e-20 ...  1.75820880e-03
         1.01517268e-03  2.69725032e-03]
       ...
       [-2.72777201e-05 -5.52598476e-04  1.75820880e-03 ... -4.87442787e-03
        -2.55311053e-03 -4.94096650e-03]
       [-6.67469497e-05 -1.13376852e-03  1.01517268e-03 ... -2.55311053e-03
        -2.56598185e-03 -3.96938994e-03]
       [-3.62691325e-05 -7.21640607e-04  2.69725032e-03 ... -4.94096650e-03
        -3.96938994e-03 -6.65773657e-03]]
    
      ...
    
      [[ 7.07882741e-03  1.66912063e-02 -2.72777201e-05 ...  1.03426213e-02
         6.77656436e-03  1.61063851e-02]
       [ 1.66912063e-02  1.33990425e-01 -3.14147778e-04 ...  1.28307748e-01
         7.49098440e-02  1.93639497e-01]
       [-2.72777201e-05 -3.14147778e-04  1.00044856e-02 ... -5.05203558e-03
        -8.52006060e-03 -6.55455060e-03]
       ...
       [ 1.03426213e-02  1.28307748e-01 -5.05203558e-03 ...  7.37378880e-01
         1.55655809e-01  6.59006666e-01]
       [ 6.77656436e-03  7.49098440e-02 -8.52006060e-03 ...  1.55655809e-01
         1.22161731e-01  2.50638198e-01]
       [ 1.61063851e-02  1.93639497e-01 -6.55455060e-03 ...  6.59006666e-01
         2.50638198e-01  9.22648331e-01]]
    
      [[ 4.64671403e-03  1.09393747e-02 -6.67469497e-05 ...  6.77656436e-03
         4.52503079e-03  1.06004080e-02]
       [ 1.09393747e-02  8.75545747e-02 -7.57438929e-04 ...  8.37531150e-02
         5.04803732e-02  1.27433073e-01]
       [-6.67469497e-05 -7.57438929e-04  6.53017747e-03 ... -3.82688183e-03
        -5.98314389e-03 -5.12386101e-03]
       ...
       [ 6.77656436e-03  8.37531150e-02 -3.82688183e-03 ...  4.80680164e-01
         1.04397439e-01  4.32931520e-01]
       [ 4.52503079e-03  5.04803732e-02 -5.98314389e-03 ...  1.04397439e-01
         8.33704088e-02  1.69134443e-01]
       [ 1.06004080e-02  1.27433073e-01 -5.12386101e-03 ...  4.32931520e-01
         1.69134443e-01  6.11230192e-01]]
    
      [[ 1.10598789e-02  2.60657172e-02 -3.62691325e-05 ...  1.61063851e-02
         1.06004080e-02  2.51616024e-02]
       [ 2.60657172e-02  2.09056153e-01 -4.17051940e-04 ...  1.99158533e-01
         1.17199073e-01  3.02285989e-01]
       [-3.62691325e-05 -4.17051940e-04  1.56084778e-02 ... -7.78106686e-03
        -1.32814228e-02 -1.01210701e-02]
       ...
       [ 1.61063851e-02  1.99158533e-01 -7.78106686e-03 ...  1.14091024e+00
         2.42430623e-01  1.02351250e+00]
       [ 1.06004080e-02  1.17199073e-01 -1.32814228e-02 ...  2.42430623e-01
         1.91322934e-01  3.92122425e-01]
       [ 2.51616024e-02  3.02285989e-01 -1.01210701e-02 ...  1.02351250e+00
         3.92122425e-01  1.44096056e+00]]]
    
    
     [[[ 1.14911426e-02  2.70372687e-02  3.16023786e-20 ...  1.66912063e-02
         1.09393747e-02  2.60657172e-02]
       [ 2.70372687e-02  2.15524108e-01  1.67989285e-18 ...  2.04720415e-01
         1.19661680e-01  3.10532067e-01]
       [ 3.16023786e-20  1.67989285e-18  1.60557730e-02 ... -7.57262516e-03
        -1.33448242e-02 -9.74746816e-03]
       ...
       [ 1.66912063e-02  2.04720415e-01 -7.57262516e-03 ...  1.16872315e+00
         2.46923793e-01  1.04751554e+00]
       [ 1.09393747e-02  1.19661680e-01 -1.33448242e-02 ...  2.46923793e-01
         1.94158375e-01  3.99191612e-01]
       [ 2.60657172e-02  3.10532067e-01 -9.74746816e-03 ...  1.04751554e+00
         3.99191612e-01  1.47362662e+00]]
    
      [[ 8.20034204e-02  2.15524108e-01  5.55697522e-19 ...  1.33990425e-01
         8.75545747e-02  2.09056153e-01]
       [ 2.15524108e-01  2.24526163e+00  2.52228894e-17 ...  2.29031800e+00
         1.31760847e+00  3.46111779e+00]
       [ 5.55697522e-19  2.52228894e-17  1.83409409e-01 ... -9.30600191e-02
        -1.61108043e-01 -1.19604516e-01]
       ...
       [ 1.33990425e-01  2.29031800e+00 -9.30600191e-02 ...  1.49372356e+01
         3.04110686e+00  1.32242040e+01]
       [ 8.75545747e-02  1.31760847e+00 -1.61108043e-01 ...  3.04110686e+00
         2.39284299e+00  4.92254063e+00]
       [ 2.09056153e-01  3.46111779e+00 -1.19604516e-01 ...  1.32242040e+01
         4.92254063e+00  1.85428598e+01]]
    
      [[ 6.88977021e-19  1.67989285e-18  1.84477971e-03 ... -3.14147778e-04
        -7.57438929e-04 -4.17051940e-04]
       [ 1.67989285e-18  2.52228894e-17  3.95275468e-02 ... -1.53878780e-02
        -2.94798027e-02 -1.99672557e-02]
       [ 1.84477971e-03  3.95275468e-02  2.88172442e-18 ...  4.70667868e-02
         2.58060369e-02  7.16660167e-02]
       ...
       [-3.14147778e-04 -1.53878780e-02  4.70667868e-02 ... -1.96175553e-01
        -9.04079287e-02 -1.93566938e-01]
       [-7.57438929e-04 -2.94798027e-02  2.58060369e-02 ... -9.04079287e-02
        -8.83774235e-02 -1.40963460e-01]
       [-4.17051940e-04 -1.99672557e-02  7.16660167e-02 ... -1.93566938e-01
        -1.40963460e-01 -2.59195437e-01]]
    
      ...
    
      [[ 7.67107360e-02  2.04720415e-01 -5.52598476e-04 ...  1.28307748e-01
         8.37531150e-02  1.99158533e-01]
       [ 2.04720415e-01  2.29031800e+00 -1.53878780e-02 ...  2.48541276e+00
         1.40932303e+00  3.67781488e+00]
       [-5.52598476e-04 -1.53878780e-02  1.93715770e-01 ... -1.26896657e-01
        -1.89009437e-01 -1.64879282e-01]
       ...
       [ 1.28307748e-01  2.48541276e+00 -1.26896657e-01 ...  1.85602186e+01
         3.61251667e+00  1.59826239e+01]
       [ 8.37531150e-02  1.40932303e+00 -1.89009437e-01 ...  3.61251667e+00
         2.79252802e+00  5.73636678e+00]
       [ 1.99158533e-01  3.67781488e+00 -1.64879282e-01 ...  1.59826239e+01
         5.73636678e+00  2.17505430e+01]]
    
      [[ 4.49946748e-02  1.19661680e-01 -1.13376852e-03 ...  7.49098440e-02
         5.04803732e-02  1.17199073e-01]
       [ 1.19661680e-01  1.31760847e+00 -2.94798027e-02 ...  1.40932303e+00
         8.85899346e-01  2.14951900e+00]
       [-1.13376852e-03 -2.94798027e-02  1.10395758e-01 ... -9.77854055e-02
        -1.30895504e-01 -1.36727518e-01]
       ...
       [ 7.49098440e-02  1.40932303e+00 -9.77854055e-02 ...  1.00883647e+01
         2.24785841e+00  9.08204180e+00]
       [ 5.04803732e-02  8.85899346e-01 -1.30895504e-01 ...  2.24785841e+00
         1.86750910e+00  3.67719192e+00]
       [ 1.17199073e-01  2.14951900e+00 -1.36727518e-01 ...  9.08204180e+00
         3.67719192e+00  1.29328073e+01]]
    
      [[ 1.16458608e-01  3.10532067e-01 -7.21640607e-04 ...  1.93639497e-01
         1.27433073e-01  3.02285989e-01]
       [ 3.10532067e-01  3.46111779e+00 -1.99672557e-02 ...  3.67781488e+00
         2.14951900e+00  5.56211761e+00]
       [-7.21640607e-04 -1.99672557e-02  2.92589671e-01 ... -1.83264236e-01
        -2.85638744e-01 -2.42279707e-01]
       ...
       [ 1.93639497e-01  3.67781488e+00 -1.83264236e-01 ...  2.65548390e+01
         5.39528232e+00  2.34166450e+01]
       [ 1.27433073e-01  2.14951900e+00 -2.85638744e-01 ...  5.39528232e+00
         4.28860356e+00  8.76246465e+00]
       [ 3.02285989e-01  5.56211761e+00 -2.42279707e-01 ...  2.34166450e+01
         8.76246465e+00  3.28996857e+01]]]
    
    
     [[[ 2.22135775e-20  3.16023786e-20  1.74187480e-04 ... -2.72777201e-05
        -6.67469497e-05 -3.62691325e-05]
       [ 3.16023786e-20  5.55697522e-19  1.84477971e-03 ... -5.52598476e-04
        -1.13376852e-03 -7.21640607e-04]
       [ 1.74187480e-04  1.84477971e-03  6.56129339e-20 ...  1.75820880e-03
         1.01517268e-03  2.69725032e-03]
       ...
       [-2.72777201e-05 -5.52598476e-04  1.75820880e-03 ... -4.87442787e-03
        -2.55311053e-03 -4.94096650e-03]
       [-6.67469497e-05 -1.13376852e-03  1.01517268e-03 ... -2.55311053e-03
        -2.56598185e-03 -3.96938994e-03]
       [-3.62691325e-05 -7.21640607e-04  2.69725032e-03 ... -4.94096650e-03
        -3.96938994e-03 -6.65773657e-03]]
    
      [[ 6.88977021e-19  1.67989285e-18  1.84477971e-03 ... -3.14147778e-04
        -7.57438929e-04 -4.17051940e-04]
       [ 1.67989285e-18  2.52228894e-17  3.95275468e-02 ... -1.53878780e-02
        -2.94798027e-02 -1.99672557e-02]
       [ 1.84477971e-03  3.95275468e-02  2.88172442e-18 ...  4.70667868e-02
         2.58060369e-02  7.16660167e-02]
       ...
       [-3.14147778e-04 -1.53878780e-02  4.70667868e-02 ... -1.96175553e-01
        -9.04079287e-02 -1.93566938e-01]
       [-7.57438929e-04 -2.94798027e-02  2.58060369e-02 ... -9.04079287e-02
        -8.83774235e-02 -1.40963460e-01]
       [-4.17051940e-04 -1.99672557e-02  7.16660167e-02 ... -1.93566938e-01
        -1.40963460e-01 -2.59195437e-01]]
    
      [[ 5.94422277e-03  1.60557730e-02  6.56129339e-20 ...  1.00044856e-02
         6.53017747e-03  1.56084778e-02]
       [ 1.60557730e-02  1.83409409e-01  2.88172442e-18 ...  1.93715770e-01
         1.10395758e-01  2.92589671e-01]
       [ 6.56129339e-20  2.88172442e-18  1.69240556e-02 ... -8.88922831e-03
        -1.51586306e-02 -1.14348051e-02]
       ...
       [ 1.00044856e-02  1.93715770e-01 -8.88922831e-03 ...  1.35884581e+00
         2.71310098e-01  1.19967636e+00]
       [ 6.53017747e-03  1.10395758e-01 -1.51586306e-02 ...  2.71310098e-01
         2.12056439e-01  4.38730120e-01]
       [ 1.56084778e-02  2.92589671e-01 -1.14348051e-02 ...  1.19967636e+00
         4.38730120e-01  1.68030330e+00]]
    
      ...
    
      [[-2.77337017e-03 -7.57262516e-03  1.75820880e-03 ... -5.05203558e-03
        -3.82688183e-03 -7.78106686e-03]
       [-7.57262516e-03 -9.30600191e-02  4.70667868e-02 ... -1.26896657e-01
        -9.77854055e-02 -1.83264236e-01]
       [ 1.75820880e-03  4.70667868e-02 -8.88922831e-03 ...  7.10845480e-02
         4.28028255e-02  1.03400107e-01]
       ...
       [-5.05203558e-03 -1.26896657e-01  7.10845480e-02 ... -1.28684924e+00
        -3.21954906e-01 -1.12856819e+00]
       [-3.82688183e-03 -9.77854055e-02  4.28028255e-02 ... -3.21954906e-01
        -2.69705680e-01 -4.98720634e-01]
       [-7.78106686e-03 -1.83264236e-01  1.03400107e-01 ... -1.12856819e+00
        -4.98720634e-01 -1.50104004e+00]]
    
      [[-4.89744402e-03 -1.33448242e-02  1.01517268e-03 ... -8.52006060e-03
        -5.98314389e-03 -1.32814228e-02]
       [-1.33448242e-02 -1.61108043e-01  2.58060369e-02 ... -1.89009437e-01
        -1.30895504e-01 -2.85638744e-01]
       [ 1.01517268e-03  2.58060369e-02 -1.51586306e-02 ...  4.28028255e-02
         3.59909013e-02  6.41591509e-02]
       ...
       [-8.52006060e-03 -1.89009437e-01  4.28028255e-02 ... -1.56940767e+00
        -3.82136493e-01 -1.42219149e+00]
       [-5.98314389e-03 -1.30895504e-01  3.59909013e-02 ... -3.82136493e-01
        -3.31307425e-01 -6.21155534e-01]
       [-1.32814228e-02 -2.85638744e-01  6.41591509e-02 ... -1.42219149e+00
        -6.21155534e-01 -2.01349114e+00]]
    
      [[-3.57050084e-03 -9.74746816e-03  2.69725032e-03 ... -6.55455060e-03
        -5.12386101e-03 -1.01210701e-02]
       [-9.74746816e-03 -1.19604516e-01  7.16660167e-02 ... -1.64879282e-01
        -1.36727518e-01 -2.42279707e-01]
       [ 2.69725032e-03  7.16660167e-02 -1.14348051e-02 ...  1.03400107e-01
         6.41591509e-02  1.56017821e-01]
       ...
       [-6.55455060e-03 -1.64879282e-01  1.03400107e-01 ... -1.63275740e+00
        -4.45122125e-01 -1.48188582e+00]
       [-5.12386101e-03 -1.36727518e-01  6.41591509e-02 ... -4.45122125e-01
        -3.91984565e-01 -7.10111521e-01]
       [-1.01210701e-02 -2.42279707e-01  1.56017821e-01 ... -1.48188582e+00
        -7.10111521e-01 -2.04561942e+00]]]
    
    
     ...
    
    
     [[[ 7.07882741e-03  1.66912063e-02 -2.72777201e-05 ...  1.03426213e-02
         6.77656436e-03  1.61063851e-02]
       [ 1.66912063e-02  1.33990425e-01 -3.14147778e-04 ...  1.28307748e-01
         7.49098440e-02  1.93639497e-01]
       [-2.72777201e-05 -3.14147778e-04  1.00044856e-02 ... -5.05203558e-03
        -8.52006060e-03 -6.55455060e-03]
       ...
       [ 1.03426213e-02  1.28307748e-01 -5.05203558e-03 ...  7.37378880e-01
         1.55655809e-01  6.59006666e-01]
       [ 6.77656436e-03  7.49098440e-02 -8.52006060e-03 ...  1.55655809e-01
         1.22161731e-01  2.50638198e-01]
       [ 1.61063851e-02  1.93639497e-01 -6.55455060e-03 ...  6.59006666e-01
         2.50638198e-01  9.22648331e-01]]
    
      [[ 7.67107360e-02  2.04720415e-01 -5.52598476e-04 ...  1.28307748e-01
         8.37531150e-02  1.99158533e-01]
       [ 2.04720415e-01  2.29031800e+00 -1.53878780e-02 ...  2.48541276e+00
         1.40932303e+00  3.67781488e+00]
       [-5.52598476e-04 -1.53878780e-02  1.93715770e-01 ... -1.26896657e-01
        -1.89009437e-01 -1.64879282e-01]
       ...
       [ 1.28307748e-01  2.48541276e+00 -1.26896657e-01 ...  1.85602186e+01
         3.61251667e+00  1.59826239e+01]
       [ 8.37531150e-02  1.40932303e+00 -1.89009437e-01 ...  3.61251667e+00
         2.79252802e+00  5.73636678e+00]
       [ 1.99158533e-01  3.67781488e+00 -1.64879282e-01 ...  1.59826239e+01
         5.73636678e+00  2.17505430e+01]]
    
      [[-2.77337017e-03 -7.57262516e-03  1.75820880e-03 ... -5.05203558e-03
        -3.82688183e-03 -7.78106686e-03]
       [-7.57262516e-03 -9.30600191e-02  4.70667868e-02 ... -1.26896657e-01
        -9.77854055e-02 -1.83264236e-01]
       [ 1.75820880e-03  4.70667868e-02 -8.88922831e-03 ...  7.10845480e-02
         4.28028255e-02  1.03400107e-01]
       ...
       [-5.05203558e-03 -1.26896657e-01  7.10845480e-02 ... -1.28684924e+00
        -3.21954906e-01 -1.12856819e+00]
       [-3.82688183e-03 -9.77854055e-02  4.28028255e-02 ... -3.21954906e-01
        -2.69705680e-01 -4.98720634e-01]
       [-7.78106686e-03 -1.83264236e-01  1.03400107e-01 ... -1.12856819e+00
        -4.98720634e-01 -1.50104004e+00]]
    
      ...
    
      [[ 4.27110568e-01  1.16872315e+00 -4.87442787e-03 ...  7.37378880e-01
         4.80680164e-01  1.14091024e+00]
       [ 1.16872315e+00  1.49372356e+01 -1.96175553e-01 ...  1.85602186e+01
         1.00883647e+01  2.65548390e+01]
       [-4.87442787e-03 -1.96175553e-01  1.35884581e+00 ... -1.28684924e+00
        -1.56940767e+00 -1.63275740e+00]
       ...
       [ 7.37378880e-01  1.85602186e+01 -1.28684924e+00 ...  3.04472750e+02
         3.79663128e+01  2.15342738e+02]
       [ 4.80680164e-01  1.00883647e+01 -1.56940767e+00 ...  3.79663128e+01
         2.63506241e+01  5.72175314e+01]
       [ 1.14091024e+00  2.65548390e+01 -1.63275740e+00 ...  2.15342738e+02
         5.72175314e+01  2.53468612e+02]]
    
      [[ 9.08275712e-02  2.46923793e-01 -2.55311053e-03 ...  1.55655809e-01
         1.04397439e-01  2.42430623e-01]
       [ 2.46923793e-01  3.04110686e+00 -9.04079287e-02 ...  3.61251667e+00
         2.24785841e+00  5.39528232e+00]
       [-2.55311053e-03 -9.04079287e-02  2.71310098e-01 ... -3.21954906e-01
        -3.82136493e-01 -4.45122125e-01]
       ...
       [ 1.55655809e-01  3.61251667e+00 -3.21954906e-01 ...  3.79663128e+01
         7.72809053e+00  3.28204336e+01]
       [ 1.04397439e-01  2.24785841e+00 -3.82136493e-01 ...  7.72809053e+00
         6.17915516e+00  1.23819197e+01]
       [ 2.42430623e-01  5.39528232e+00 -4.45122125e-01 ...  3.28204336e+01
         1.23819197e+01  4.51089103e+01]]
    
      [[ 3.83671838e-01  1.04751554e+00 -4.94096650e-03 ...  6.59006666e-01
         4.32931520e-01  1.02351250e+00]
       [ 1.04751554e+00  1.32242040e+01 -1.93566938e-01 ...  1.59826239e+01
         9.08204180e+00  2.34166450e+01]
       [-4.94096650e-03 -1.93566938e-01  1.19967636e+00 ... -1.12856819e+00
        -1.42219149e+00 -1.48188582e+00]
       ...
       [ 6.59006666e-01  1.59826239e+01 -1.12856819e+00 ...  2.15342738e+02
         3.28204336e+01  1.67865912e+02]
       [ 4.32931520e-01  9.08204180e+00 -1.42219149e+00 ...  3.28204336e+01
         2.40911149e+01  5.13305925e+01]
       [ 1.02351250e+00  2.34166450e+01 -1.48188582e+00 ...  1.67865912e+02
         5.13305925e+01  2.15342738e+02]]]
    
    
     [[[ 4.64671403e-03  1.09393747e-02 -6.67469497e-05 ...  6.77656436e-03
         4.52503079e-03  1.06004080e-02]
       [ 1.09393747e-02  8.75545747e-02 -7.57438929e-04 ...  8.37531150e-02
         5.04803732e-02  1.27433073e-01]
       [-6.67469497e-05 -7.57438929e-04  6.53017747e-03 ... -3.82688183e-03
        -5.98314389e-03 -5.12386101e-03]
       ...
       [ 6.77656436e-03  8.37531150e-02 -3.82688183e-03 ...  4.80680164e-01
         1.04397439e-01  4.32931520e-01]
       [ 4.52503079e-03  5.04803732e-02 -5.98314389e-03 ...  1.04397439e-01
         8.33704088e-02  1.69134443e-01]
       [ 1.06004080e-02  1.27433073e-01 -5.12386101e-03 ...  4.32931520e-01
         1.69134443e-01  6.11230192e-01]]
    
      [[ 4.49946748e-02  1.19661680e-01 -1.13376852e-03 ...  7.49098440e-02
         5.04803732e-02  1.17199073e-01]
       [ 1.19661680e-01  1.31760847e+00 -2.94798027e-02 ...  1.40932303e+00
         8.85899346e-01  2.14951900e+00]
       [-1.13376852e-03 -2.94798027e-02  1.10395758e-01 ... -9.77854055e-02
        -1.30895504e-01 -1.36727518e-01]
       ...
       [ 7.49098440e-02  1.40932303e+00 -9.77854055e-02 ...  1.00883647e+01
         2.24785841e+00  9.08204180e+00]
       [ 5.04803732e-02  8.85899346e-01 -1.30895504e-01 ...  2.24785841e+00
         1.86750910e+00  3.67719192e+00]
       [ 1.17199073e-01  2.14951900e+00 -1.36727518e-01 ...  9.08204180e+00
         3.67719192e+00  1.29328073e+01]]
    
      [[-4.89744402e-03 -1.33448242e-02  1.01517268e-03 ... -8.52006060e-03
        -5.98314389e-03 -1.32814228e-02]
       [-1.33448242e-02 -1.61108043e-01  2.58060369e-02 ... -1.89009437e-01
        -1.30895504e-01 -2.85638744e-01]
       [ 1.01517268e-03  2.58060369e-02 -1.51586306e-02 ...  4.28028255e-02
         3.59909013e-02  6.41591509e-02]
       ...
       [-8.52006060e-03 -1.89009437e-01  4.28028255e-02 ... -1.56940767e+00
        -3.82136493e-01 -1.42219149e+00]
       [-5.98314389e-03 -1.30895504e-01  3.59909013e-02 ... -3.82136493e-01
        -3.31307425e-01 -6.21155534e-01]
       [-1.32814228e-02 -2.85638744e-01  6.41591509e-02 ... -1.42219149e+00
        -6.21155534e-01 -2.01349114e+00]]
    
      ...
    
      [[ 9.08275712e-02  2.46923793e-01 -2.55311053e-03 ...  1.55655809e-01
         1.04397439e-01  2.42430623e-01]
       [ 2.46923793e-01  3.04110686e+00 -9.04079287e-02 ...  3.61251667e+00
         2.24785841e+00  5.39528232e+00]
       [-2.55311053e-03 -9.04079287e-02  2.71310098e-01 ... -3.21954906e-01
        -3.82136493e-01 -4.45122125e-01]
       ...
       [ 1.55655809e-01  3.61251667e+00 -3.21954906e-01 ...  3.79663128e+01
         7.72809053e+00  3.28204336e+01]
       [ 1.04397439e-01  2.24785841e+00 -3.82136493e-01 ...  7.72809053e+00
         6.17915516e+00  1.23819197e+01]
       [ 2.42430623e-01  5.39528232e+00 -4.45122125e-01 ...  3.28204336e+01
         1.23819197e+01  4.51089103e+01]]
    
      [[ 7.15552870e-02  1.94158375e-01 -2.56598185e-03 ...  1.22161731e-01
         8.33704088e-02  1.91322934e-01]
       [ 1.94158375e-01  2.39284299e+00 -8.83774235e-02 ...  2.79252802e+00
         1.86750910e+00  4.28860356e+00]
       [-2.56598185e-03 -8.83774235e-02  2.12056439e-01 ... -2.69705680e-01
        -3.31307425e-01 -3.91984565e-01]
       ...
       [ 1.22161731e-01  2.79252802e+00 -2.69705680e-01 ...  2.63506241e+01
         6.17915516e+00  2.40911149e+01]
       [ 8.33704088e-02  1.86750910e+00 -3.31307425e-01 ...  6.17915516e+00
         5.34563196e+00  1.03009602e+01]
       [ 1.91322934e-01  4.28860356e+00 -3.91984565e-01 ...  2.40911149e+01
         1.03009602e+01  3.51186004e+01]]
    
      [[ 1.46801848e-01  3.99191612e-01 -3.96938994e-03 ...  2.50638198e-01
         1.69134443e-01  3.92122425e-01]
       [ 3.99191612e-01  4.92254063e+00 -1.40963460e-01 ...  5.73636678e+00
         3.67719192e+00  8.76246465e+00]
       [-3.96938994e-03 -1.40963460e-01  4.38730120e-01 ... -4.98720634e-01
        -6.21155534e-01 -7.10111521e-01]
       ...
       [ 2.50638198e-01  5.73636678e+00 -4.98720634e-01 ...  5.72175314e+01
         1.23819197e+01  5.13305925e+01]
       [ 1.69134443e-01  3.67719192e+00 -6.21155534e-01 ...  1.23819197e+01
         1.03009602e+01  2.05466755e+01]
       [ 3.92122425e-01  8.76246465e+00 -7.10111521e-01 ...  5.13305925e+01
         2.05466755e+01  7.42144278e+01]]]
    
    
     [[[ 1.10598789e-02  2.60657172e-02 -3.62691325e-05 ...  1.61063851e-02
         1.06004080e-02  2.51616024e-02]
       [ 2.60657172e-02  2.09056153e-01 -4.17051940e-04 ...  1.99158533e-01
         1.17199073e-01  3.02285989e-01]
       [-3.62691325e-05 -4.17051940e-04  1.56084778e-02 ... -7.78106686e-03
        -1.32814228e-02 -1.01210701e-02]
       ...
       [ 1.61063851e-02  1.99158533e-01 -7.78106686e-03 ...  1.14091024e+00
         2.42430623e-01  1.02351250e+00]
       [ 1.06004080e-02  1.17199073e-01 -1.32814228e-02 ...  2.42430623e-01
         1.91322934e-01  3.92122425e-01]
       [ 2.51616024e-02  3.02285989e-01 -1.01210701e-02 ...  1.02351250e+00
         3.92122425e-01  1.44096056e+00]]
    
      [[ 1.16458608e-01  3.10532067e-01 -7.21640607e-04 ...  1.93639497e-01
         1.27433073e-01  3.02285989e-01]
       [ 3.10532067e-01  3.46111779e+00 -1.99672557e-02 ...  3.67781488e+00
         2.14951900e+00  5.56211761e+00]
       [-7.21640607e-04 -1.99672557e-02  2.92589671e-01 ... -1.83264236e-01
        -2.85638744e-01 -2.42279707e-01]
       ...
       [ 1.93639497e-01  3.67781488e+00 -1.83264236e-01 ...  2.65548390e+01
         5.39528232e+00  2.34166450e+01]
       [ 1.27433073e-01  2.14951900e+00 -2.85638744e-01 ...  5.39528232e+00
         4.28860356e+00  8.76246465e+00]
       [ 3.02285989e-01  5.56211761e+00 -2.42279707e-01 ...  2.34166450e+01
         8.76246465e+00  3.28996857e+01]]
    
      [[-3.57050084e-03 -9.74746816e-03  2.69725032e-03 ... -6.55455060e-03
        -5.12386101e-03 -1.01210701e-02]
       [-9.74746816e-03 -1.19604516e-01  7.16660167e-02 ... -1.64879282e-01
        -1.36727518e-01 -2.42279707e-01]
       [ 2.69725032e-03  7.16660167e-02 -1.14348051e-02 ...  1.03400107e-01
         6.41591509e-02  1.56017821e-01]
       ...
       [-6.55455060e-03 -1.64879282e-01  1.03400107e-01 ... -1.63275740e+00
        -4.45122125e-01 -1.48188582e+00]
       [-5.12386101e-03 -1.36727518e-01  6.41591509e-02 ... -4.45122125e-01
        -3.91984565e-01 -7.10111521e-01]
       [-1.01210701e-02 -2.42279707e-01  1.56017821e-01 ... -1.48188582e+00
        -7.10111521e-01 -2.04561942e+00]]
    
      ...
    
      [[ 3.83671838e-01  1.04751554e+00 -4.94096650e-03 ...  6.59006666e-01
         4.32931520e-01  1.02351250e+00]
       [ 1.04751554e+00  1.32242040e+01 -1.93566938e-01 ...  1.59826239e+01
         9.08204180e+00  2.34166450e+01]
       [-4.94096650e-03 -1.93566938e-01  1.19967636e+00 ... -1.12856819e+00
        -1.42219149e+00 -1.48188582e+00]
       ...
       [ 6.59006666e-01  1.59826239e+01 -1.12856819e+00 ...  2.15342738e+02
         3.28204336e+01  1.67865912e+02]
       [ 4.32931520e-01  9.08204180e+00 -1.42219149e+00 ...  3.28204336e+01
         2.40911149e+01  5.13305925e+01]
       [ 1.02351250e+00  2.34166450e+01 -1.48188582e+00 ...  1.67865912e+02
         5.13305925e+01  2.15342738e+02]]
    
      [[ 1.46801848e-01  3.99191612e-01 -3.96938994e-03 ...  2.50638198e-01
         1.69134443e-01  3.92122425e-01]
       [ 3.99191612e-01  4.92254063e+00 -1.40963460e-01 ...  5.73636678e+00
         3.67719192e+00  8.76246465e+00]
       [-3.96938994e-03 -1.40963460e-01  4.38730120e-01 ... -4.98720634e-01
        -6.21155534e-01 -7.10111521e-01]
       ...
       [ 2.50638198e-01  5.73636678e+00 -4.98720634e-01 ...  5.72175314e+01
         1.23819197e+01  5.13305925e+01]
       [ 1.69134443e-01  3.67719192e+00 -6.21155534e-01 ...  1.23819197e+01
         1.03009602e+01  2.05466755e+01]
       [ 3.92122425e-01  8.76246465e+00 -7.10111521e-01 ...  5.13305925e+01
         2.05466755e+01  7.42144278e+01]]
    
      [[ 5.40067365e-01  1.47362662e+00 -6.65773657e-03 ...  9.22648331e-01
         6.11230192e-01  1.44096056e+00]
       [ 1.47362662e+00  1.85428598e+01 -2.59195437e-01 ...  2.17505430e+01
         1.29328073e+01  3.28996857e+01]
       [-6.65773657e-03 -2.59195437e-01  1.68030330e+00 ... -1.50104004e+00
        -2.01349114e+00 -2.04561942e+00]
       ...
       [ 9.22648331e-01  2.17505430e+01 -1.50104004e+00 ...  2.53468612e+02
         4.51089103e+01  2.15342738e+02]
       [ 6.11230192e-01  1.29328073e+01 -2.01349114e+00 ...  4.51089103e+01
         3.51186004e+01  7.42144278e+01]
       [ 1.44096056e+00  3.28996857e+01 -2.04561942e+00 ...  2.15342738e+02
         7.42144278e+01  3.04472750e+02]]]]


## Now we can build Fock Matrix and Solve the generalized eigenvalue problem


```python
from scipy.linalg import eigh
nelectrons = 1+1+8 #H H O for H2O
norb = nelectrons // 2
maxiter = 100
C = np.eye(n) #coeffcient matrix
prev = float('inf')
tol = 1e-7
for i in range(maxiter):
    #start building Fock matrix
    JK = np.zeros((n,n))
    F = np.zeros((n,n)) 
    for j in range(n):
        for k in range(n):
            for l in range(n):
                for m in range(n):
                    for o in range(norb):
                        JK[j][k] += (C[l][o]*C[m][o] * 
                                    (2*coulombRepulsionTensor[j][k][l][m]-
                                     coulombRepulsionTensor[j][m][k][l]
                                    ))
    F = Tmatrix + coulombAttractionMatrix + JK #add kinetic integral and coulomb attraction integeral
    S = Smatrix #overlap matrix
    #now solve the FC=SCe (Hartree Fock Roothaan-Hall Equation)
    energy, C = eigh(F, S)
    EHF = 0
    for i in range(norb):
        EHF += 2*energy[i] - 1/2*JK[i][i]
    if abs(EHF - prev) < tol:
        print('SCF Converged')
        break
    delta = EHF - prev
    prev = EHF
    print('EHF:'+str(EHF)+" "+'prev:'+str(prev)+' '+'delta:'+str(delta))
```

    EHF:-89.0867859100989 prev:-89.0867859100989 delta:-inf
    EHF:-67.12761944636007 prev:-67.12761944636007 delta:21.959166463738825
    EHF:-67.6195916882164 prev:-67.6195916882164 delta:-0.49197224185633104
    EHF:-67.66855653433281 prev:-67.66855653433281 delta:-0.04896484611640517
    EHF:-67.67630737732246 prev:-67.67630737732246 delta:-0.007750842989651119
    EHF:-67.67778511136507 prev:-67.67778511136507 delta:-0.0014777340426093133
    EHF:-67.67810581339319 prev:-67.67810581339319 delta:-0.00032070202811951276
    EHF:-67.67817983523364 prev:-67.67817983523364 delta:-7.402184044735804e-05
    EHF:-67.67819739475206 prev:-67.67819739475206 delta:-1.7559518425969145e-05
    EHF:-67.67820161464026 prev:-67.67820161464026 delta:-4.219888197098953e-06
    EHF:-67.67820263601112 prev:-67.67820263601112 delta:-1.0213708634410068e-06
    EHF:-67.67820288436762 prev:-67.67820288436762 delta:-2.4835649981014285e-07
    SCF Converged


## Now let's visualize the MO
1. First Build MO List from Converged C
2. Normalize the MO
3. Visualize MO (HOMO)

We can see for water, the homo is mainly Oxygen's P orbital.


```python
molist = []
for i in range(n):
    coeffs = C[:,i]
    newmo = Mo(aolist,coeffs)
    molist.append(newmo)
plot_implicit(molist[4], 0.01,[-5,5],0,90) #choose HOMO
```


![png](output_26_0.png)



```python

```
