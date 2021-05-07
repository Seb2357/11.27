# -------------------------------------------------------------------------
#
# PYTHON for DUMMIES 20-21
# Problème 8
#
# Script de test
#  Vincent Legat
#
# Largement inspiré du programme de Nicolas Roisin :-)
# Ou les méthodes numériques pour obtenir la solution du projet P2 !
#
# -------------------------------------------------------------------------
#


from numpy import *
from scipy.spatial import Delaunay

# ------------------------------------------------------------------------------------
#
# Intégration du flux magnétique
#
#    Xmagnet,Ymagnet : tableaux numpy des coordonnées x,z des sommets  de l'aimant 
#    Zmagnet : scalaire contenant la hauteur du l'aimant
#    Xcoil,Ycoil : tableaux numpy des coordonnées x,z des sommets de la bobine
#    triangles : tableau contenant les indices des 3 sommets de chaque élément
#    Xshift : tableau numpy contenant les translation de l'aimant sur une période
#    mu0 : perméabilité du vide
#    mu  : valeur absolue de la composante z du momemt magnétique du l'aimant = [0,0,-mu]
#   
#  La fonction renvoie un vecteur phi contenant le flux du champs magnétique intercepté
#  par une spire exprimé en [T cm2]
#

def area(x,y) :
    
    a = sqrt((x[0]-x[1])**2+(y[0]-y[1])**2)
    b = sqrt((x[0]-x[2])**2+(y[0]-y[2])**2)
    c = sqrt((x[1]-x[2])**2+(y[1]-y[2])**2)
    
    p = (a+b+c)/2
    
    area = sqrt(p*(p-a)*(p-b)*(p-c))
    
    return area

def X_I(i, Xmagnet, Ymagnet, triangles) :
    
    x_i = mean(Xmagnet[triangles[i,:]])
    y_i = mean(Ymagnet[triangles[i,:]])
    
    return x_i,y_i
    
def r(a,b,c,Xshift): 
    
    """ 
    a est le point d'arrivée du vecteur et b est le point de départ (qui est mobile)
    """
    
    r = sqrt((a[0]-(b[0]-Xshift))**2 + (a[1]-b[1])**2 + c**2)
    
    return r


def magnetComputeInduction(Xmagnet,Ymagnet,Zmagnet,Xcoil,Ycoil,triangles,Xshift,mu0,mu) :
  
    m = len(Xshift)
    phi = zeros(m)
    n = len(triangles)
    total_area = 0
    
    for i in range(n):
        total_area += area(Xmagnet[triangles[i,:]], Ymagnet[triangles[i,:]])
        
    Area_coil = zeros(n)
    Area_magnet = zeros(n)
    for i in range(n):
        Area_coil[i] = area(Xcoil[triangles[i,:]], Ycoil[triangles[i,:]])
        Area_magnet[i] = area(Xmagnet[triangles[i,:]], Ymagnet[triangles[i,:]])
        
    X_I_coil = zeros((2,n))
    X_I_magnet = zeros((2,n))
    for i in range(n):
        X_I_coil[0][i], X_I_coil[1][i] = X_I(i, Xcoil, Ycoil, triangles)
        X_I_magnet[0][i], X_I_magnet[1][i] = X_I(i, Xmagnet, Ymagnet, triangles)
        
    R = zeros((m,n,n))
    for i in range(m):
        for j in range(n):
            for k in range(n):
                R[i][j][k] = r(X_I_coil.T[j],X_I_magnet.T[k],Zmagnet,Xshift[i])
    
    for i in range(m):              #i-ème intervalle de temps
        for j in range(n):          #j-ème élément de la bobine
            for k in range(n):      #k-ème élément de l'aimant
                phi[i] += Area_coil[j] * ((mu0/(4*pi*R[i][j][k]**3))*((3*-mu*(Area_magnet[k] / total_area)*(Zmagnet/abs((R[i][j][k])))*(Zmagnet/abs((R[i][j][k]))))-(-mu*(Area_magnet[k] / total_area))))
        print(" Iteration %2d  : shift = %6.3f [cm] : phi = %.8f" % (i,Xshift[i],phi[i]))        
        
    return phi
 
# ------------------------------------------------------------------------------------ 
#
# Script de test
#
#
# -0- Paramètres matériels
#
# ------------------------------------------------------------------------------------


mu0     = 4e-7*pi*1e-2     # permeabilité du vide en [H/cm] 
Rmagnet = 1.25             # rayon de l'aimant [cm]
Hmagnet = 0.6              # épaisseur de l'aimant [cm]
Zmagnet = 0.5              # position verticale de l'aimant en [cm]
Br      = 1.4              # magnetisation residuelle du Néodyme fer bore (NdFeB) en [T] ou [kg/(A s)]
mu      = Rmagnet**2*Hmagnet*pi*Br / mu0    
                           # moment magnétique de l'aimant [A cm2]
Rcoil   = 1                # rayon de la bobine [cm]
nSpires = 200



# ------------------------------------------------------------------------------------
#
# -1- Construction d'un maillage de triangles pour un cercle de rayon unitaire
#
# ------------------------------------------------------------------------------------


nR      = 6
nTheta  = 6
nNode   = 1 + sum(arange(1,nR))*nTheta
R     = zeros(nNode)
Theta = zeros(nNode)

index = 1; dR = 1.0/(nR-1)
for i in range(1,nR):
    dTheta = 2*pi/(i*nTheta)
    for j in range(0,i*nTheta):
        R[index]     = i*dR
        Theta[index] = j*dTheta; index += 1

X       = R*cos(Theta)
Y       = R*sin(Theta)

triangles = Delaunay(stack((X,Y),1)).simplices
nElem = len(triangles)

print(" Number of triangles : %d " % nElem)
print(" Number of nodes     : %d " % nNode)


# ------------------------------------------------------------------------------------
#
# -2- Calcul du flux et de la tension induite dans la bobine
#
# ------------------------------------------------------------------------------------

m       = 41
Xstart  = -5                        # [cm]
Xstop   =  5                        # [cm]
Xshift  = linspace(Xstart,Xstop,m)
Tstart  = 0                         # [s]
Tstop   = 0.5                       # [s]
T,delta = linspace(Tstart,Tstop,m,retstep=True)

Xmagnet = Rmagnet*R*cos(Theta)
Ymagnet = Rmagnet*R*sin(Theta) 
Xcoil   = Rcoil*R*cos(Theta)
Ycoil   = Rcoil*R*sin(Theta) 
    
phi     = magnetComputeInduction(Xmagnet,Ymagnet,Zmagnet,Xcoil,Ycoil,triangles,
                                                                   Xshift,mu0,mu)  
phi     = phi * nSpires    
voltage = - diff(phi) / (delta*10)

# ------------------------------------------------------------------------------------
#
# -3- Quelques jolis plots et animation
#
# ------------------------------------------------------------------------------------

import matplotlib.pyplot as plt
import matplotlib.animation as animation
plt.rcParams['toolbar'] = 'None'

def frame(i):
  plt.clf()

  n = 50
  X,Z = meshgrid(linspace(-2,2,n),linspace(-2,2,n))
  Y  = zeros_like(X)
  Bx = zeros(shape(X))
  Bz = zeros(shape(X))

  for iElem in range(nElem):
    Xp = X - Xdipole[iElem] - Xshift[i]
    Yp = Y - Ydipole[iElem]
    Zp = Z - Zmagnet
    r     = sqrt(Xp*Xp + Yp*Yp + Zp*Zp)
    coeff = -(mu0*mu) / (4*pi*r**5)
    Bx   += coeff * (3*Zp*Xp)
    Bz   += coeff * (3*Zp*Zp - r*r)
  plt.streamplot(X,Z,Bx,Bz, density=1.4, linewidth=None, color='blue')

  x = array([-Rmagnet,Rmagnet,Rmagnet,-Rmagnet,-Rmagnet]) + Xshift[i]
  y = array([0,0,Hmagnet,Hmagnet,0])+Zmagnet-Hmagnet/2.0
  plt.fill(x,y,facecolor='blue',alpha=1)

  x = [-Rcoil,Rcoil]
  y = [0,0]
  plt.plot(x,y,"-r",linewidth=4)
  
  plt.xlim((-2,2)); plt.ylim((-2,2))
  plt.title('Electromagnetic Field')   

# ------------------------------------------------------------------------------------

fig=plt.figure("Maillage de l'aimant")
plt.plot(Xmagnet,Ymagnet,'or')
plt.triplot(Xmagnet,Ymagnet,triangles,'-k')
Xdipole = mean(Xmagnet[triangles[:,:]],axis=1)  
Ydipole = mean(Ymagnet[triangles[:,:]],axis=1)  
plt.plot(Xdipole,Ydipole,'ob')  

plt.axis("equal")
plt.axis("off")

plt.figure("Flux et tension induite sur une période")
plt.plot(T,phi,'-r')
plt.plot(T[1:],voltage,'-b')
plt.text(0.01,-100,"$N\phi(t)$ [T cm$^2$]",color='red',fontsize=12)
plt.text(0.3,100,r"$-N\dfrac{\partial \phi}{\partial t}(t)$ [mV]",color='blue',fontsize=12)
plt.text(0.4,-210,"time [s]",color='black',fontsize=12)

plt.figure("Un joli plot pour le coordinateur :-)",figsize=(10, 10))
frame(20)

movie = animation.FuncAnimation(plt.figure("Claude's project",figsize=(10,10)),frame,41,interval=20,repeat=False)
plt.show()









