import numpy as np
import matplotlib.pylab as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import mpl_toolkits.mplot3d.art3d as a3

################################################################################
################################################################################
def get_collisions(infile,verbose=False):

    collisions = []

    not_at_end = True
    collision_count = 0
    new_collision = True
    while ( not_at_end ):

        ############################################################################
        # Read in one collision
        ############################################################################
        line = infile.readline()

        if collision_count%1000==0 and verbose:
            print "collision count: ",collision_count

        if line=="":
            not_at_end = False

        if line.find("Event")>=0:
            new_collision = True

        if new_collision==True:

            # Read in the pion info for this collision.
            pions = []
            line = infile.readline()
            npions = int(line)
            for i in xrange(npions):
                line = infile.readline()
                vals = line.split()
                e = float(vals[0])
                px = float(vals[1])
                py = float(vals[2])
                pz = float(vals[3])
                q = int(vals[4])
                sigpi = float(vals[5])
                sigka = float(vals[6])
                likpi = float(vals[7])
                likka = float(vals[8])
                nphopi = int(vals[9])
                nphoka = int(vals[10])
                depthmu = float(vals[11])
                cluster_energy = float(vals[12])
                pions.append([e,px,py,pz,q,sigpi,sigka,likpi,likka,nphopi,nphoka,depthmu,cluster_energy])

            # Read in the kaon info for this collision.
            kaons = []
            line = infile.readline()
            nkaons = int(line)
            for i in xrange(nkaons):
                line = infile.readline()
                vals = line.split()
                e = float(vals[0])
                px = float(vals[1])
                py = float(vals[2])
                pz = float(vals[3])
                q = int(vals[4])
                sigpi = float(vals[5])
                sigka = float(vals[6])
                likpi = float(vals[7])
                likka = float(vals[8])
                nphopi = int(vals[9])
                nphoka = int(vals[10])
                depthmu = float(vals[11])
                cluster_energy = float(vals[12])
                kaons.append([e,px,py,pz,q,sigpi,sigka,likpi,likka,nphopi,nphoka,depthmu,cluster_energy])

            # Read in the muon info for this collision.
            muons = []
            line = infile.readline()
            nmuons = int(line)
            for i in xrange(nmuons):
                line = infile.readline()
                vals = line.split()
                e = float(vals[0])
                px = float(vals[1])
                py = float(vals[2])
                pz = float(vals[3])
                q = int(vals[4])
                sigpi = float(vals[5])
                sigka = float(vals[6])
                likpi = float(vals[7])
                likka = float(vals[8])
                nphopi = int(vals[9])
                nphoka = int(vals[10])
                depthmu = float(vals[11])
                cluster_energy = float(vals[12])
                muons.append([e,px,py,pz,q,sigpi,sigka,likpi,likka,nphopi,nphoka,depthmu,cluster_energy])

            # Read in the electron info for this collision.
            electrons = []
            line = infile.readline()
            nelectrons = int(line)
            for i in xrange(nelectrons):
                line = infile.readline()
                vals = line.split()
                e = float(vals[0])
                px = float(vals[1])
                py = float(vals[2])
                pz = float(vals[3])
                q = int(vals[4])
                sigpi = float(vals[5])
                sigka = float(vals[6])
                likpi = float(vals[7])
                likka = float(vals[8])
                nphopi = int(vals[9])
                nphoka = int(vals[10])
                depthmu = float(vals[11])
                cluster_energy = float(vals[12])
                electrons.append([e,px,py,pz,q,sigpi,sigka,likpi,likka,nphopi,nphoka,depthmu,cluster_energy])


            # Read in the photon info for this collision.
            photons = []
            line = infile.readline()
            nphotons = int(line)
            for i in xrange(nphotons):
                line = infile.readline()
                vals = line.split()
                e = float(vals[0])
                px = float(vals[1])
                py = float(vals[2])
                pz = float(vals[3])
                photons.append([e,px,py,pz])


            # Read in the information about the missing transverse energy (MET) in the collision.
            # This is really the x and y direction for the missing momentum.
            #line = infile.readline()
            #vals = line.split()
            #met_px = float(vals[0])
            #met_py = float(vals[1])

            new_collision = False
            collision_count += 1

            collisions.append([pions,kaons,muons,electrons,photons])

    return collisions

################################################################################
################################################################################
def draw_jet(origin=(0,0),angle=90,length=0.5,opening_angle=20,ntracks=5,show_tracks=False):

    lines = []
    patches = []

    # Edges of cone
    width_at_top = length*np.deg2rad(opening_angle)
    for side in [-1,1]:
        theta0 = np.deg2rad(angle+(side*opening_angle/2.0)) 
        x1 = length*np.cos(theta0)
        y1 = length*np.sin(theta0)
        #print x1,y1
        line = mlines.Line2D((origin[0],x1), (origin[1],y1), lw=2., alpha=0.4,color='red',markeredgecolor='red')
        lines.append(line)

    # End of cone
    arad = np.deg2rad(angle)
    center = (origin[0]+np.cos(arad)*length,origin[1]+np.sin(arad)*length)
    #print center
    p = mpatches.Ellipse(center, width_at_top+0.01, width_at_top/2.0,facecolor='red',alpha=0.4,edgecolor='gray',angle=abs(angle+90))
    patches.append(p)

    return patches,lines


    

################################################################################
################################################################################
def draw_jets(origins=[(0,0)],angles=[90],lengths=[0.5],opening_angles=[20],ntrackss=[5],show_trackss=[False]):

    alllines = []
    allpatches = []

    # Edges of cone
    for origin,angle,length,opening_angle,ntracks,show_tracks in zip(origins,angles,lengths,opening_angles,ntrackss,show_trackss):
        patches,lines = draw_jet(origin=origin,angle=angle,length=length,opening_angle=opening_angle,ntracks=ntracks,show_tracks=show_tracks)
        allpatches += patches
        alllines += lines


    return allpatches,alllines


    
################################################################################
################################################################################
def draw_line3D(origin=[(0,0,0)],pmom=[(1,1,1)],color='red',ls='-',lw=2.0):

    lines = []

    #print pmom
    for o,p in zip(origin,pmom):
        #x1 = p[0]
        #y1 = p[1]
        #z1 = p[2]
        x1 = p[2]
        y1 = p[0]
        z1 = p[1]
        #print x1,y1,z1
        line = a3.Line3D((o[0],x1),(o[1],y1),(o[0],z1), lw=lw, ls=ls, alpha=0.9,color=color,markeredgecolor=color)
        lines.append(line)

    return lines


################################################################################
################################################################################
def draw_beams():

    lines = draw_line3D(origin=[(0,0,-0.1),(0,0,0.1)],pmom=[(0,0,-1.0),(0,0,1.0)],color='green',lw=1)

    return lines

################################################################################
################################################################################
def draw_jet3D(origin=[(0,0,0)],pmom=[(1,1,1)]):

    neworg = origin.copy()
    newmom = pmom.copy()

    offset = [[0.05,0.05,0.05],
              [0.05,0.05,-0.05],
              [0.05,-0.05,0.05],
              [0.05,-0.05,-0.05],
              [-0.05,0.05,0.05],
              [-0.05,0.05,-0.05],
              [-0.05,-0.05,0.05],
              [-0.05,-0.05,-0.05],
            ]

    offset = np.array(offset)
    offset *= 50

    for p in pmom:
        for o in offset:
            #print p.copy(),o
            pnew = p.copy() + o
            #print pnew
            newmom = np.vstack((newmom,pnew))
            neworg = np.vstack((neworg,(0,0,0)))

    lines = draw_line3D(origin=neworg,pmom=newmom,color='orange',lw=1)

    return lines

################################################################################
################################################################################
def draw_pion3D(origin=[(0,0,0)],pmom=[(1,1,1)]):

    lines = draw_line3D(origin=origin,pmom=pmom,color='red',lw=5)

    return lines

################################################################################
def draw_kaon3D(origin=[(0,0,0)],pmom=[(1,1,1)]):

    lines = draw_line3D(origin=origin,pmom=pmom,color='orange',lw=5)

    return lines

################################################################################
################################################################################
def draw_muon3D(origin=[(0,0,0)],pmom=[(1,1,1)]):

    lines = draw_line3D(origin=origin,pmom=pmom,color='blue',lw=5)

    return lines

################################################################################
################################################################################
def draw_electron3D(origin=[(0,0,0)],pmom=[(1,1,1)]):

    lines = draw_line3D(origin=origin,pmom=pmom,color='green',lw=2)

    return lines


################################################################################
################################################################################
def draw_photon3D(origin=[(0,0,0)],pmom=[(1,1,1)]):

    lines = draw_line3D(origin=origin,pmom=pmom,color='gray',ls='-',lw=4)

    return lines

################################################################################
################################################################################
def display_collision3D(collision):

    pions,kaons,muons,electrons,photons = collision

    lines = draw_beams()

    pmom = np.array(pions).transpose()[1:4].transpose()
    origin = np.zeros((len(pions),3))
    lines += draw_pion3D(origin=origin,pmom=pmom)

    pmom = np.array(kaons).transpose()[1:4].transpose()
    origin = np.zeros((len(kaons),3))
    lines += draw_kaon3D(origin=origin,pmom=pmom)

    pmom = np.array(muons).transpose()[1:4].transpose()
    origin = np.zeros((len(muons),3))
    lines += draw_muon3D(origin=origin,pmom=pmom)

    pmom = np.array(electrons).transpose()[1:4].transpose()
    origin = np.zeros((len(electrons),3))
    lines += draw_electron3D(origin=origin,pmom=pmom)

    pmom = np.array(photons).transpose()[1:4].transpose()
    origin = np.zeros((len(photons),3))
    lines += draw_photon3D(origin=origin,pmom=pmom)

    fig = plt.figure(figsize=(7,5),dpi=100)
    ax = fig.add_subplot(1,1,1)
    ax = fig.gca(projection='3d')
    plt.subplots_adjust(top=0.98,bottom=0.02,right=0.98,left=0.02)

    for l in lines:
        ax.add_line(l)

    ax.set_xlim(-1,1)
    ax.set_ylim(-1,1)
    ax.set_zlim(-1,1)

    #return lines,fig,ax

