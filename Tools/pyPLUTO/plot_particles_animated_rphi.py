from matplotlib import animation
from pylab import *
import matplotlib.colors as colors
import pyPLUTO.pload as pp # importing the pyPLUTO pload module.
import pyPLUTO.ploadparticles as pr # importing the pyPLUTO ploadparticles module.
from matplotlib.animation import FuncAnimation
def plot_particles_animated_rphi(ntot, w_dir, UNIT_DENSITY, UNIT_LENGTH, UNIT_VELOCITY, xmin, xmax, datatype, file_name = 'particles_rphi.gif', out_dir = ""):
    f1 = plt.figure(figsize=[10,8])
    plt.rcParams["figure.dpi"] = 200
    plt.rcParams['axes.linewidth'] = 0.1
    P = pr.ploadparticles(0, w_dir, datatype=datatype,ptype='CR') # Loading particle data : particles.00ns_ch00.flt

    PVmag = np.sqrt(P.vx1**2 + P.vx2**2 + P.vx3**2) # estimating the velocity magnitude
    maxU = 0
    if(len(PVmag > 0)):
        max(PVmag)
    index = 0
    for i in range(ntot+1):
        P = pr.ploadparticles(i, w_dir, datatype=datatype, ptype='CR')
        PVmag = np.sqrt(P.vx1 ** 2 + P.vx2 ** 2 + P.vx3 ** 2)  # estimating the velocity magnitude
        if(len(PVmag) > 0):
            u = max(PVmag)
            if(u > maxU):
                maxU = u
                index = i

    #index = 27
    P = pr.ploadparticles(index, w_dir, datatype=datatype, ptype='CR')
    PVmag = np.sqrt(P.vx1 ** 2 + P.vx2 ** 2 + P.vx3 ** 2)

    particles = np.zeros((len(P.x1), 2))
    for i in range(len(particles)):
        particles[i][0] = P.x1[i]*sin(P.x2[i])*cos(P.x3[i])*UNIT_LENGTH
        particles[i][1] = P.x1[i]*sin(P.x2[i])*sin(P.x3[i])*UNIT_LENGTH
        
    ### background magnetic field
    D = pp.pload(ntot, varNames=['Bx1', 'Bx2', 'Bx3'], w_dir=w_dir, datatype=datatype)  # Load fluid data.
    ndim = len((D.Bx1.shape))

    minB = 0
    maxB = 0
    nx = 0
    ny = 0
    
    nx = D.Bx1.shape[0]

    if (ndim == 1):
        ny = 10
    else :
        ny = D.Bx1.shape[1]

    B = np.zeros([ny, nx])
    
    if (ndim == 1):
        Bz = D.Bx3.T[:] * np.sqrt(4 * np.pi * UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY)
        By = D.Bx2.T[:] * np.sqrt(4 * np.pi * UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY)
        Bx = D.Bx1.T[:] * np.sqrt(4 * np.pi * UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY)
        for i in range(ny):
            B[i] = np.sqrt(np.square(Bx) + np.square(By) + np.square(Bz))
    if (ndim == 2):
        Bz = D.Bx3.T[:, :] * np.sqrt(4 * np.pi * UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY)
        By = D.Bx2.T[:, :] * np.sqrt(4 * np.pi * UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY)
        Bx = D.Bx1.T[:, :] * np.sqrt(4 * np.pi * UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY)
        B = np.sqrt(np.square(Bx) + np.square(By) + np.square(Bz))
    if (ndim == 3):
        zpoint = math.floor(D.Bx1.T.shape[0] / 2)
        Bz = D.Bx3.T[zpoint, :, :] * np.sqrt(4 * np.pi * UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY)
        By = D.Bx2.T[zpoint, :, :] * np.sqrt(4 * np.pi * UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY)
        Bx = D.Bx1.T[zpoint, :, :] * np.sqrt(4 * np.pi * UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY)
        B = np.sqrt(np.square(Bx) + np.square(By) + np.square(Bz))
    #np.flip(B, 0)

    minB = np.amin(B)
    maxB = np.amax(B)
    
    xmin = D.x1.min() * UNIT_LENGTH
    xmax = D.x1.max() * UNIT_LENGTH

    for i in range(ntot + 1):
        D = pp.pload(i, varNames=['Bx1', 'Bx2', 'Bx3'], w_dir=w_dir, datatype=datatype)  # Load fluid data.
        if (ndim == 1):
            Bz = D.Bx3.T[:] * np.sqrt(4 * np.pi * UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY)
            By = D.Bx2.T[:] * np.sqrt(4 * np.pi * UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY)
            Bx = D.Bx1.T[:] * np.sqrt(4 * np.pi * UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY)
            for i in range(ny):
                B[i] = np.sqrt(np.square(Bx) + np.square(By) + np.square(Bz))
        if (ndim == 2):
            Bz = D.Bx3[:, :] * np.sqrt(4 * np.pi * UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY)
            By = D.Bx2[:, :] * np.sqrt(4 * np.pi * UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY)
            Bx = D.Bx1[:, :] * np.sqrt(4 * np.pi * UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY)
            B = np.sqrt(np.square(Bx) + np.square(By) + np.square(Bz))
        if (ndim == 3):
            zpoint = math.floor(D.Bx1.T.shape[0] / 2)
            Bz = D.Bx3[zpoint, :, :] * np.sqrt(4 * np.pi * UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY)
            By = D.Bx2[zpoint, :, :] * np.sqrt(4 * np.pi * UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY)
            Bx = D.Bx1[zpoint, :, :] * np.sqrt(4 * np.pi * UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY)
            B = np.sqrt(np.square(Bx) + np.square(By) + np.square(Bz))
        if (np.amin(B) < minB):
            minB = np.amin(B)
        if (np.amax(B) > maxB):
            maxB = np.amax(B)

    print("maxB = ", maxB)
    print("minB = ", minB)


    def update(frame_number):
        f1.clear()
        ax = f1.add_subplot(projection="polar")
        cax1 = f1.add_axes([0.91, 0.12, 0.03, 0.75])
        cax2 = f1.add_axes([0.125, 0.92, 0.75, 0.03])
        P = pr.ploadparticles(frame_number, w_dir, datatype=datatype,
                              ptype='CR')  # Loading particle data : particles.00ns_ch00.flt

        PVmag = np.sqrt(P.vx1 ** 2 + P.vx2 ** 2 + P.vx3 ** 2)  # estimating the velocity magnitude

        particles = np.zeros((len(P.x1), 2))
        for i in range(len(particles)):
            #particles[i][0] = P.x1[i]*sin(P.x2[i])*cos(P.x3[i])*UNIT_LENGTH
            #particles[i][1] = P.x1[i]*sin(P.x2[i])*sin(P.x3[i])*UNIT_LENGTH
            particles[i][1] = P.x1[i]*UNIT_LENGTH
            particles[i][0] = P.x3[i]

        #ax.set_xlim([-xmax, xmax])
        #ax.set_ylim([-xmax, xmax])
        ax.set_title('Number of particles = ' + str(len(particles)))
        
        B = np.zeros([ny, nx])
        
        D = pp.pload(frame_number, varNames=['Bx1', 'Bx2', 'Bx3'], w_dir=w_dir, datatype=datatype)  # Load fluid data.
        if (ndim == 1):
            Bz = D.Bx3.T[:] * np.sqrt(4 * np.pi * UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY)
            By = D.Bx2.T[:] * np.sqrt(4 * np.pi * UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY)
            Bx = D.Bx1.T[:] * np.sqrt(4 * np.pi * UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY)
            for i in range(ny):
                B[i] = np.sqrt(np.square(Bx) + np.square(By) + np.square(Bz))
        if (ndim == 2):
            Bz = D.Bx3[:, :] * np.sqrt(4 * np.pi * UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY)
            By = D.Bx2[:, :] * np.sqrt(4 * np.pi * UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY)
            Bx = D.Bx1[:, :] * np.sqrt(4 * np.pi * UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY)
            B = np.sqrt(np.square(Bx) + np.square(By) + np.square(Bz))
        if (ndim == 3):
            zpoint = math.floor(D.Bx1.shape[2] / 2)
            Bz = D.Bx3[zpoint, :, :] * np.sqrt(4 * np.pi * UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY)
            By = D.Bx2[zpoint, :, :] * np.sqrt(4 * np.pi * UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY)
            Bx = D.Bx1[zpoint, :, :] * np.sqrt(4 * np.pi * UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY)
            B = np.sqrt(np.square(Bx) + np.square(By) + np.square(Bz))

        rad = np.linspace(0, xmax, nx)
        azm = np.linspace(0, 2 * np.pi, ny)
        r, th = np.meshgrid(rad, azm)
        im2 = ax.pcolormesh(th, r, B)
        plt.colorbar(im2, cax=cax2, orientation='horizontal', norm=colors.Normalize(vmin=minB, vmax=maxB))

        #im1 = ax.scatter(particles[:,0], particles[:,1], s=10, c=PVmag, cmap=colors.Normalize(vmin=0, vmax=maxU))  # scatter plot
        im1 = ax.scatter(particles[:,0], particles[:,1], s=10, c=PVmag, cmap=plt.get_cmap('hot'), vmin = 0, vmax = maxU)  # scatter plot
        plt.colorbar(im1, cax=cax1)  # vertical colorbar for particle data.
            

            
        #plt.subplot(projection="polar")

        
        
        
        return im1

    anim = FuncAnimation(f1, update, interval=10, frames = ntot+1)

    #plt.show()

    f = out_dir + file_name
    writergif = animation.PillowWriter(fps=4)
    anim.save(f, writer=writergif)
    plt.close()
