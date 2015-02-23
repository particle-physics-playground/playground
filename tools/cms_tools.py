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

            # Read in the jet info for this collision.
            jets = []
            line = infile.readline()
            njets = int(line)
            for i in xrange(njets):
                line = infile.readline()
                vals = line.split()
                e = float(vals[0])
                px = float(vals[1])
                py = float(vals[2])
                pz = float(vals[3])
                bquark_jet_tag = float(vals[4])
                jets.append([e,px,py,pz,bquark_jet_tag])

            # Read in the muon info for this collision.
            muons = []
            line = infile.readline()
            nmuons = int(line)
            num_mu=0
            for i in xrange(nmuons):
                line = infile.readline()
                vals = line.split()
                e = float(vals[0])
                px = float(vals[1])
                py = float(vals[2])
                pz = float(vals[3])
                charge = int(vals[4])
                muons.append([e,px,py,pz,charge])
                num_mu+=1
                

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
                charge = int(vals[4])
                electrons.append([e,px,py,pz,charge])

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
            line = infile.readline()
            vals = line.split()
            met_px = float(vals[0])
            met_py = float(vals[1])

            new_collision = False
            collision_count += 1

            collisions.append([jets,muons,electrons,photons,[met_px,met_py]])

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

    lines = draw_line3D(origin=[(0,0,-0.1),(0,0,0.1)],pmom=[(0,0,-200.0),(0,0,200.0)],color='red',lw=1)

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

    lines = draw_line3D(origin=origin,pmom=pmom,color='gray',ls='--',lw=7)

    return lines

################################################################################
################################################################################
def display_collision3D(collision):

    jets,muons,electrons,photons,met = collision

    lines = draw_beams()

    pmom = np.array(jets).transpose()[1:4].transpose()
    origin = np.zeros((len(jets),3))
    lines += draw_jet3D(origin=origin,pmom=pmom)

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

    ax.set_xlim(-200,200)
    ax.set_ylim(-200,200)
    ax.set_zlim(-200,200)

    #return lines,fig,ax

