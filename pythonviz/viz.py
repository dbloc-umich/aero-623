import numpy as np
import matplotlib.pyplot as plt
from readgri import readgri
import sys


#-----------------------------------------------------------
def readU(fname):
    return np.loadtxt(fname)

    
#-----------------------------------------------------------
def cp(P):
    g=1.4
    gmi = 1.4-1
    P0=1/g
    Pout=P0*0.7
    Mout2 = 2/gmi*((P0/Pout)**(gmi/g)-1)
    qout = 0.5*g*Pout*Mout2
    cp = (P-Pout)/qout

    return cp


#-----------------------------------------------------------
def getEdgeP(Mesh,U):
    valsUpper=[]
    valsLower=[]
    V = Mesh['V']; B_up = Mesh['B_up']; B_low=Mesh['B_low']
    U = getField(U,'mach')
    Pup=U[B_up[:,2],5]
    Plow=U[B_low[:,2],5]
    xup = (V[B_up[:,0],0]+V[B_up[:,1],0])/2 # x location of first node
    xlow = (V[B_low[:,0],0]+V[B_low[:,1],0])/2 # x location of first node

    return Pup, Plow, xup, xlow

#-----------------------------------------------------------
def plotCoefficients(Mesh,U):
    Pup, Plow, xup, xlow = getEdgeP(Mesh,U)
    cp_up = cp(Pup)
    cp_low = cp(Plow)
    fig, ax = plt.subplots(figsize=(8, 6))

    ax.plot(xup,  cp_up,  color='blue',  label='Upper Surface')
    ax.plot(xlow, cp_low, color='red',   label='Lower Surface')

    ax.set_title('Pressure Coefficient Distribution', fontsize=14)
    ax.set_xlabel('x/c', fontsize=12)
    ax.set_ylabel('$C_p$', fontsize=12)
    ax.legend(fontsize=12)
    ax.invert_yaxis()  # Cp plots conventionally have -Cp pointing up
    ax.grid(True, linestyle='--', alpha=0.5)

    plt.tight_layout()
    plt.savefig('pressureCoeff.png', dpi=150)
    plt.show()

    return 0
#-----------------------------------------------------------
def getField(U, field):
    r, ru, rv, rE = [U[:,i] for i in range(4)]
    g = 1.4
    V = np.sqrt(ru**2 + rv**2)/r
    p = (g-1.)*(rE-0.5*r*V**2)
    c = np.sqrt(g*p/r)
    M = V/c
    S = p/r**g
    U = np.column_stack([U, V, p, c, M, S])
    # print(U.shape)
    q = field.lower()
    if (q == 'mach'):
        return U
    else:
        return []

#-----------------------------------------------------------
def getnodestate(Mesh, U):
    V = Mesh['V']; E = Mesh['E']
    Nv, Ne = V.shape[0], E.shape[0]
    UN = np.zeros(Nv); count = np.zeros(Nv)
    for e in range(Ne):
        for i in range(3):
            n = E[e,i];
            UN[n] += U[e]; count[n] += 1
    UN /= count
    return UN

#-----------------------------------------------------------
def plotmesh(Mesh, fname):
    V = Mesh['V']; E = Mesh['E']; BE = Mesh['BE']
    f = plt.figure(figsize=(12,12))
    plt.triplot(V[:,0], V[:,1], E, 'k-')
    for i in range(BE.shape[0]):
        plt.plot(V[BE[i,0:2],0],V[BE[i,0:2],1], '-', linewidth=1, color='black')
    dosave = not not fname
    plt.axis('equal'); plt.axis('off')
    #plt.axis([-0.5, 1.5,-1, 1])
    plt.tick_params(axis='both', labelsize=12)
    f.tight_layout();
    if (dosave): plt.savefig(fname)
    else: plt.show(block=True);
    plt.close(f)
    
#-----------------------------------------------------------
def plotstate(Mesh, U, p, field, frange, fname):
    V = Mesh['V']; E = Mesh['E']; BE = Mesh['BE']
    f = plt.figure(figsize=(12,12))
    U_mod = getField(U, field)
    # print(U_mod.shape)
    F=U_mod[:,7] # plotting Mach number
    if (p == 0):
        plt.tripcolor(V[:,0], V[:,1], triangles=E, facecolors=F, shading='flat')
    else:
        vc = np.linspace(frange[0], frange[1], 21) if (len(frange) > 0) else 20
        plt.tricontourf(V[:,0], V[:,1], E, getnodestate(Mesh,F), vc)
    for i in range(BE.shape[0]):
        plt.plot(V[BE[i,0:2],0],V[BE[i,0:2],1], '-', linewidth=2, color='black')
    dosave = not not fname
    plt.axis('equal'); plt.axis('off')
    if (len(frange)>0): plt.clim(frange[0], frange[1])
    plt.set_cmap('jet')
    cbar=plt.colorbar(orientation='horizontal', pad=-.12, fraction=.045)
    cbar.ax.tick_params(labelsize=16)
    #plt.axis([-0.5, 1.5,-1, 1])
    plt.tick_params(axis='both', labelsize=12)
    f.tight_layout();
    if (dosave): plt.savefig(fname, bbox_inches='tight',pad_inches=-.2)
    else: plt.show(block=True);
    plt.close(f)


#-----------------------------------------------------------
def main():
    if (len(sys.argv) < 2):
        print('Pass at least one argument: Mesh')
    else:
        # fMesh, fU, p, field, fname
        Mesh = readgri(sys.argv[1])
        U = [] if (len(sys.argv) <= 2) else readU(sys.argv[2])
        p = 0 if (len(sys.argv) <= 3) else int(sys.argv[3])
        field = 'Mach' if (len(sys.argv) <= 4) else sys.argv[4]
        fname = [] if (len(sys.argv) <= 5) else sys.argv[5]
        if len(U) > 0:
            frange = []
            if (field.lower() == 'mach'):
                frange = [0.0, 0.5]
            plotstate(Mesh, U, p, field, frange, fname)
            plotCoefficients(Mesh,U)
        else:
            plotmesh(Mesh, 'mesh.png')
    
if __name__ == "__main__":
    main()

