#!/usr/bin/python
import numpy, math, time, argparse, sys, pickle, struct, message_pb2
import writing_binary, writing_pdb
import random as random
from constants import EMPTY, BINDER, LAMIN, BSITE_R, BSITE_L, REGDNA, DIST

#accepted move vectors
#MOVES = numpy.array([[0,0,1], [0,1,0], [1,0,0], [0,0,-1], [0,-1,0], [-1,0,0], [1,0,1], [0,1,1], [1,1,0], [-1,0,-1], [0,-1,-1], [-1,-1,0], [1,0,-1], [0,1,-1], [1,-1,0], [-1,0,1], [0,-1,1], [-1,1,0], [-1,1,-1], [-1,-1,-1], [1,1,1], [1,-1,1], [1,-1,-1], [-1,-1,1], [-1,1,1], [1,1,-1]])
MOVES = numpy.array([[0,0,1], [0,1,0], [1,0,0], [0,0,-1], [0,-1,0], [-1,0,0], [1,0,1], [0,1,1], [1,1,0], [-1,0,-1], [0,-1,-1], [-1,-1,0], [1,0,-1], [0,1,-1], [1,-1,0], [-1,0,1], [0,-1,1], [-1,1,0]])
#accepted matching positions of binding sites
#BMOVES = numpy.array([[1,0,0], [-1,0,0], [0,1,0], [0,-1,0], [0,0,1], [0,0,-1]])
BMOVES = numpy.array([[1,0,0], [-1,0,0], [0,1,0], [0,-1,0], [0,0,1], [0,0,-1], [1,1,0], [-1,-1,0], [-1,1,0], [1,-1,0], [1,0,1], [-1,0,1], [1,0,-1], [-1, 0, -1], [0,1,-1], [0,-1,-1], [0,-1,1], [0,1,1]])
#TERMOVES = numpy.array([[-2, -2, -2], [-2, -2, -1], [-2, -2, 0], [-2, -2, 1], [-2, -2, 2], [-2, -1, -2], [-2, -1, -1], [-2, -1, 0], [-2, -1, 1], [-2, -1, 2], [-2, 0, -2], [-2, 0, -1], [-2, 0, 0], [-2, 0, 1], [-2, 0, 2], [-2, 1, -2], [-2, 1, -1], [-2, 1, 0], [-2, 1, 1], [-2, 1, 2], [-2, 2, -2], [-2, 2, -1], [-2, 2, 0], [-2, 2, 1], [-2, 2, 2], [-1, -2, -2], [-1, -2, -1], [-1, -2, 0], [-1, -2, 1], [-1, -2, 2], [-1, -1, -2], [-1, -1, -1], [-1, -1, 0], [-1, -1, 1], [-1, -1, 2], [-1, 0, -2], [-1, 0, -1], [-1, 0, 0], [-1, 0, 1], [-1, 0, 2], [-1, 1, -2], [-1, 1, -1], [-1, 1, 0], [-1, 1, 1], [-1, 1, 2], [-1, 2, -2], [-1, 2, -1], [-1, 2, 0], [-1, 2, 1], [-1, 2, 2], [0, -2, -2], [0, -2, -1], [0, -2, 0], [0, -2, 1], [0, -2, 2], [0, -1, -2], [0, -1, -1], [0, -1, 0], [0, -1, 1], [0, -1, 2], [0, 0, -2], [0, 0, -1], [0, 0, 1], [0, 0, 2], [0, 1, -2], [0, 1, -1], [0, 1, 0], [0, 1, 1], [0, 1, 2], [0, 2, -2], [0, 2, -1], [0, 2, 0], [0, 2, 1], [0, 2, 2], [1, -2, -2], [1, -2, -1], [1, -2, 0], [1, -2, 1], [1, -2, 2], [1, -1, -2], [1, -1, -1], [1, -1, 0], [1, -1, 1], [1, -1, 2], [1, 0, -2], [1, 0, -1], [1, 0, 0], [1, 0, 1], [1, 0, 2], [1, 1, -2], [1, 1, -1], [1, 1, 0], [1, 1, 1], [1, 1, 2], [1, 2, -2], [1, 2, -1], [1, 2, 0], [1, 2, 1], [1, 2, 2], [2, -2, -2], [2, -2, -1], [2, -2, 0], [2, -2, 1], [2, -2, 2], [2, -1, -2], [2, -1, -1], [2, -1, 0], [2, -1, 1], [2, -1, 2], [2, 0, -2], [2, 0, -1], [2, 0, 0], [2, 0, 1], [2, 0, 2], [2, 1, -2], [2, 1, -1], [2, 1, 0], [2, 1, 1], [2, 1, 2], [2, 2, -2], [2, 2, -1], [2, 2, 0], [2, 2, 1], [2, 2, 2]])
TERMOVES = numpy.array([[1,0,0], [-1,0,0], [0,1,0], [0,-1,0], [0,0,1], [0,0,-1], [1,1,0], [-1,-1,0], [-1,1,0], [1,-1,0], [1,0,1], [-1,0,1], [1,0,-1], [-1, 0, -1], [0,1,-1], [0,-1,-1], [0,-1,1], [0,1,1]])

# radius of the nucleus
R = 20
# 2 x radius + a fringe, because lamin barrier has to be hermetic
BOUND = 2 * R + 2
GOOD_NEIGH = 3

def pars_inp():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", dest="Out_str", default='out', required=False,
                        help="The prefix of the output filenames with lamins and with MC trajectory.")
    parser.add_argument('--binary', action='store_true', help="Whether to use binary output format.")
    parser.add_argument("-i", dest="In_str", default = '', required=False,
                        help="An input state (pickled file).")
    parser.add_argument("-s", type=int, dest="Steps", default=100000, required=False,
                        help="Number of steps for simulation (default 100000).")
    parser.add_argument("-l",  dest="Ch_lenght", default=512, required=False,
                        help="Length of the chains, separated by comma (two chains: 234, 456) (default one chain with 512 bins lenght).")
    parser.add_argument("-r", dest="Revers", default=1058, required=False,
                        help="Length of the chains, where tail should be near the center (2R, 3R), separated by comma (two chains: 234, 456) (default one chain with 512 bins lenght).")
    parser.add_argument("-b", type=int, dest="Binders", default=256, required=False,
                        help="Number of binders (default 256).")
    parser.add_argument("-a","--laminBSitesFile", required=False, default="",
                        help="A file with positions of binding sites for lamins.")
    parser.add_argument("-n", "--regularBSitesFile", required=False, default = "",
                        help="A file with positions of binding sites of floating binders.")
    parser.add_argument("--simulation_name", required=False, default="Example simulation",
                        help="Simulation name.")
    parser.add_argument("--simulation_description", required=False, default="Example description",
                        help="Simulation description.")
    parser.add_argument("--fpkf", type=int, default=10, required=False,
                        help="""Number of frames per keyframe (default 10),
                                equal to the number of delta frames + 1.
                                A frame is outputed after each fpkf - 1 moves.""")

    args = parser.parse_args()
    if (args.laminBSitesFile == "" or args.regularBSitesFile == "") and args.In_str == "":
        print "Error, if you do not provide the -i option, then you have to provide the --laminBSitesFile and --regularBSitesFile options"
        sys.exit(1)
    return args	

def init_dist_matrix(max_d = GOOD_NEIGH + 1):
    dist_matrix = numpy.zeros((max_d, max_d, max_d), dtype = numpy.float32)
    for i in range(max_d):
        for j in range(max_d):
            for k in range(max_d):
                dist_matrix[i][j][k] = math.sqrt(i**2 + j**2 + k**2)
    return dist_matrix
DIST_MATRIX = init_dist_matrix()

def initialize_import(f):
    list_ob = pickle.load(open(f))
    print len(list_ob)
    ch = list_ob[0]
    b = list_ob[1]
    a = list_ob[2]
    state = list_ob[3]
    ch_n = list_ob[4]
    re = list_ob[5]
    return ch, b, a, state, ch_n, re

def dist_from_mi(x, y, z, mi):
    """Return the Euclidean distance of (x, y, z) from (mi, mi, mi)."""
    return math.sqrt((x - mi)**2 + (y - mi)**2 + (z - mi)**2)
    
def getStateWithLamins(bound, f, binary):
    """Return a 3D grid of states with lamins positioned on a sphere and write the states into a file.""" 
    state = numpy.zeros((bound, bound, bound), dtype=numpy.int)
    MIDDLE = bound / 2
    binders_list = []
    for x in range(BOUND):
        for y in range(BOUND):
            for z in range(BOUND):
                border = abs(dist_from_mi(x, y, z, MIDDLE) - MIDDLE + 1)
                if border <= 2:
                    state[x, y, z] = LAMIN
                    if border == 1:
                        binders_list.append([x, y , z])
    if binary:
        lam_name = f.split('.bin')[0] + '_lamin.bin'
        writing_binary.write_lamins_binary(lam_name, binders_list)
    else:
        lam_name = f.split('.pdb')[0] + '_lamin.pdb'
        writing_pdb.write_lamins_pdb(lam_name, binders_list)
    return state

def getTerritories(angles):
    """Return a state with territories for each chain. 

    The keyword angles should be in the format [(0,20), (20,70), (70,200), (200,360)]."""
    stateT = numpy.zeros((BOUND, BOUND, BOUND), dtype = numpy.int)
    
    radius_cut = round(MIDDLE/5.0)
    start = 0.
    radius_cut_list = []
    for ra in range(5):
        start = start+radius_cut
        radius_cut_list.append(start)
            
    for katy in angles: 
        kat0 =  katy[0]
        kat1 = katy[1]

        for x in range(BOUND):
            for y in range(BOUND):
                for z in range(BOUND):
                    border = dist_from_mi(x-MIDDLE, y-MIDDLE, z-MIDDLE, 0)
                    if border <= MIDDLE and border!= 0:
                        kat_tet = math.atan2((y-MIDDLE),(x-MIDDLE))
                        # If kat_tet is between angles kat0 and kat1
                        if (kat0 < 180 and kat1 < 180 and math.radians(kat0) <= kat_tet < math.radians(kat1))\
                                or (kat0 < 180 and kat1 > 180\
                                    and ((math.radians(kat0) <= kat_tet <= math.pi)\
                                         or (-math.pi < kat_tet < (-math.pi+math.radians(kat1-180)))))\
                                or (kat0 > 180 and kat1 > 180\
                                    and (-math.pi+math.radians(kat0-180) < kat_tet < -math.pi+math.radians(kat1 - 180))):
                            bor = dist_from_mi(x - MIDDLE, y - MIDDLE, MIDDLE - MIDDLE, 0)
                            if 0. <= bor < radius_cut_list[0]: 
                                stateT[x, y, z] = angles.index(katy) + 1
                            elif radius_cut_list[0] <= bor < radius_cut_list[1]: 
                                stateT[x, y, z] = str(angles.index(katy)+1) + '1'
                            elif radius_cut_list[1] <= bor < radius_cut_list[2]: 
                                stateT[x, y, z] = str(angles.index(katy)+1) + '2'
                            elif radius_cut_list[2] <= bor < radius_cut_list[3]: 
                                stateT[x, y, z] = str(angles.index(katy)+1) + '3'
                            elif radius_cut_list[3] <= bor < MIDDLE: 
                                stateT[x, y, z] = str(angles.index(katy)+1) + '4'
    return stateT

def dist(r1,r2):
    return DIST_MATRIX[tuple(abs(r1 - r2))]
    
#def distance (p1, p2):
#   return numpy.sqrt(numpy.sum((p1-p2)**2))

def cross(cha, po1, po2):
    pier = numpy.where((cha == po1).all(axis=1))
    dru = numpy.where((cha == po2).all(axis=1))
    #print "pir-dru ", po1, po2, len(pier), len(dru), pier, dru, pier[0][0]
    if len(pier) == 1 and len(dru) == 1:
        if pier[0][0] == dru[0][0]+1 or pier[0][0] == dru[0][0]-1:
            #print 'TRUE cross'
            return True
        else:
            #print 'FALSE cross'
            return False
    else: print "Two atoms of the chain have the same coordinates. Please check ", po1, po2, pier[0], dru[0]
    
def intersect(new_p, next, sta, ch): #coordinates of two next polimer points
    #print dist(next, new_p), distance(next, new_p)
    #print "INTER"
    if  dist(next, new_p) == 1: 
        return False
    elif dist(next, new_p) > 1:
        differ = new_p - next
        if differ[0] == 0:
            pos1 = numpy.array(([next[0], next[1], new_p[2]]))
            pos2 = numpy.array(([next[0], new_p[1], next[2]]))
            if sta[tuple(pos1)] not in [EMPTY, LAMIN, BINDER] and sta[tuple(pos2)] not in [EMPTY, LAMIN, BINDER]:
                return cross(ch, pos1, pos2)
            else: return False         
        elif differ[1] == 0:
            pos1 = numpy.array(([next[0], next[1], new_p[2]]))
            pos2 = numpy.array(([new_p[0], next[1], next[2]]))
            if sta[tuple(pos1)] not in [EMPTY, LAMIN, BINDER] and sta[tuple(pos2)] not in [EMPTY, LAMIN, BINDER]:
                return cross(ch, pos1, pos2)
            else: return False
        elif differ[2] == 0:
            pos1 = numpy.array(([next[0], new_p[1], next[2]]))
            pos2 = numpy.array(([new_p[0], next[1], next[2]]))
            if sta[tuple(pos1)] not in [EMPTY, LAMIN, BINDER] and sta[tuple(pos2)] not in [EMPTY, LAMIN, BINDER]:
                return cross(ch, pos1, pos2)
            else: return False
        else: print new_p, next, "The distance between these two positions are not = sqrt(2)"

def no_collisions(x, state):
    """Return True if the state in x is empty and False otherwise."""
    #print "COLI", x, state.shape, state[tuple(x)]
    if state[tuple(x)] != EMPTY:
        return False
    else:
        return True
        
def get_angles(at_nr):
    """Return a list of pairs of angles of beginnings and ends of consecutive chains

    at_nr is a list of lengths of consecutive chains
    """
    stos = [round(float(a)/min(at_nr)) for a in at_nr]
    ang_st = 360./sum(stos)
    ang = [round(s * ang_st) for s in stos]
    mi = min(ang)
    if mi < 15:
        rest = 15-mi
        ma = max(ang)
        ang[ang.index(mi)] = 15
        ang[ang.index(ma)] = ang[ang.index(ma)] - rest
    an0 = ang[0]
    angles = [(0, ang[0])]
    for an in ang[1:-1]:
        an1 = an0+an
        angles.append((an0, an1))
        an0 = an1
    angles.append((an0, 360))
    return angles
    
MIDDLE = BOUND / 2
def initialize_random(n, m, fa, rever, regularBSitesFile, laminBSitesFile, binary, bound = BOUND): 
    """Initialize the state randomly.
    
    n - a list with chains atoms numbers, m - number of binders, fa - file name of lamins
    """
    def fibonacci_sphere(samples=1, randomize=True): 
        """Distribute points on a sphere from web."""
        rnd = 1.
        if randomize:
            rnd = random.random() * samples
        points = []
        offset = 2./samples
        increment = math.pi * (3. - math.sqrt(5.));
        for i in range(samples):
            y = round(((i * offset) - 1) + (offset / 2))
            r = math.sqrt(9 - pow(y,2)) # R ~ 3
            #print "R", r
            phi = ((i + rnd) % samples) * increment
            y = y + MIDDLE
            x = round(math.cos(phi) * r) + MIDDLE
            z = round(math.sin(phi) * r) + MIDDLE
            points.append([x,y,z])
            #print points
        return points

    def get_site_type_list(fpath, length_list):
        positions = []
        chrom = -1
        #print "lent", length_list
        for length in length_list:
            pos = [0] * length
            positions.append(pos)
        for l in open(fpath):
                if "chr" in l:
                    chrom +=1
                    continue
                positions[chrom][int(l) -1] = 1
        return positions
        
    def get_site_type(i, regular_bsites, lamin_bsites): # BSITE_R interacts with binders whereas BSITE_L interacts both with lamins and binders
        if lamin_bsites[i] == 1:
            return BSITE_L
        elif regular_bsites[i] == 1:
            return BSITE_R
        else:
            return REGDNA
    
    def check_gyr(atoms, nr_a, chai):
        r = 0.0
        for at in range(atoms-nr_a, atoms+1):
            for ato in range(at, atoms+1):
                r = r + math.sqrt(numpy.sum((chai[at] - chai[ato])**2))
                #print at, ato, r, chai[at], chai[ato], (chai[at] - chai[ato])**2,numpy.sum((chai[at] - chai[ato])**2) 
        ch_gyr = 2*r/float(nr_a)
        #print "NUMERY", atoms-nr_a, atoms+1, nr_a, ch_gyr
        return ch_gyr
        
    def rand_next(cu, st, ch, terri, nun, ti, at =1, ii=1):
        mov = random.choice(MOVES)
        tries = 0
        #print "dfgf"
        ch_cop = numpy.copy(ch)
        ch_cop[at]=cu + mov
        #while (tries < 100) and (not (no_collisions(tuple(cu + mov), st)) or intersect(cu, cu+mov, st, ch)  or terri[tuple(cu+mov)] != 0 or check_gyr(at, ii, ch_cop)>nun/5.):
        while (tries < 100) and (not (no_collisions(tuple(cu + mov), st)) or intersect(cu, cu + mov, st, ch)  or terri[tuple(cu + mov)] != ti):
        #while tries < 100 and (not (no_collisions(tuple(cu + mov), st)) or intersect(cu, cu+mov, st, ch)):
            #print "TRTRT"
            mov = random.choice(MOVES)
            tries += 1
            try:
                terri[tuple(cu+mov)]
            except IndexError: 
                print "Error"
                mov = random.choice(MOVES)
            ch_cop[at]=cu + mov
        if tries == 100: 
            return mov, True
        return mov, False
        
    def fill_one_terr(one_t, pos, t):
        for tm in TERMOVES:
            try:
                if one_t[tuple(tm + pos)] == t:
                    one_t[tuple(tm + pos)] = 0
            except IndexError: pass
        return one_t
    
    chain = numpy.zeros((sum(n), 3), dtype = numpy.int)
    binders = numpy.zeros((m, 3), dtype = numpy.int)
    state = getStateWithLamins(bound, fa, binary)
    ter_ang = get_angles(n)
    territor = getTerritories(ter_ang)
    
    attached_to_lamins = []
    for x in range(bound):
        for y in range(bound):
            for z in range(bound):
                dist_m = dist_from_mi(x, y, z, MIDDLE)
                if dist_m <= 3:
                    territor[x, y, z] = 0

    regular_bsites = get_site_type_list(regularBSitesFile, n)
    lamin_bsites   = get_site_type_list(laminBSitesFile, n)
    #print regular_bsites      
    
    at_nr = -1
    
    #points_on_sph = fibonacci_sphere(60)
    #print get_site_type(0, regular_bsites[0], lamin_bsites[0])
    for nu, re, la in zip(n, regular_bsites, lamin_bsites):
         
        cur0 = [bound / 2] * 3
        te = n.index(nu)+1
        if cur0 == [MIDDLE,MIDDLE,MIDDLE]:
            ter_point = numpy.where(territor == te)
            po_dis_min = 20.0 # big enough to find smaller dis - closest to the middle
            for px, py, pz in zip(ter_point[0], ter_point[1], ter_point[2]):
                po_dis = dist_from_mi(px, py, pz, MIDDLE)
                if po_dis_min > po_dis:
                    po_dis_min = po_dis
                    cur0 = [px,py,pz]
        
        at_nr += 1
        chain[at_nr] = numpy.array(cur0)
        print chain[0], chain[at_nr]
        if nu not in rever:
            state[tuple(chain[at_nr])] = get_site_type(0, re, la)
        else:
            state[tuple(chain[at_nr])] = get_site_type(-1, re, la)
        cur = chain[at_nr]
        print at_nr, chain[at_nr], nu
        
        ter_zero = numpy.where(territor == te)
        ter_one = numpy.where(territor == int(str(te) + '1'))
        ter_two = numpy.where(territor == int(str(te) + '2'))
        ter_three = numpy.where(territor == int(str(te) + '3'))
        ter_four = numpy.where(territor == int(str(te) + '4'))
        state_lam = numpy.where(state == LAMIN)
        try_chain = 0 
        succ = False  
        
        if nu not in rever:
            while try_chain < 100 and not succ:
                print try_chain, succ
                for i in range(1, nu):
                    chain_l_dev = round(nu/125.)
                    at_nr += 1
                    if chain_l_dev <= i < chain_l_dev * 10. and len(numpy.where(territor == int(str(te) + '1'))[0]) != 0:
                        print "ZMIENIAM", str(te) + '1'
                        open_part = numpy.where(territor == int(str(te) + '1'))
                        territor[open_part]= te
                    elif chain_l_dev * 10. <= i < chain_l_dev * 30. and len(numpy.where(territor == int(str(te) + '2'))[0]) != 0:
                        print "ZMIENIAM", str(te) + '2'
                        open_part = numpy.where(territor == int(str(te) + '2'))
                        territor[open_part]= te
                    elif chain_l_dev * 30. <= i < chain_l_dev * 70. and len(numpy.where(territor == int(str(te) + '3'))[0]) != 0:
                        print "ZMIENIAM", str(te) + '3'
                        open_part = numpy.where(territor == int(str(te) + '3'))
                        territor[open_part]= te
                    elif chain_l_dev * 70. <= i < chain_l_dev * 125. and len(numpy.where(territor == int(str(te) + '4'))[0]) != 0:
                        print "ZMIENIAM", str(te) + '4'
                        open_part = numpy.where(territor == int(str(te) + '4'))
                        territor[open_part]= te 
            
                    #print "AT_nr!!!", at_nr, i,  cur, type(cur)
                    mo, tang = rand_next(cur, state, chain, territor, nu, te, at_nr, i)
                    if tang:
                        try_chain +=1
                        assert try_chain != 100, "Chain is not possible to initialize"
                        at_nr = at_nr - i
                        old_state0 = state[tuple(chain[at_nr])]
                        territor[ter_zero]= te
                        state[ter_zero] = 0.
                        territor[ter_one]= int(str(te) + '1')
                        state[ter_one] = 0.
                        territor[ter_two]= int(str(te) + '2')
                        state[ter_two] = 0.
                        territor[ter_three]= int(str(te) + '3')
                        state[ter_three] = 0.
                        territor[ter_four]= int(str(te) + '4')
                        state[ter_four] = 0.
                        state[tuple(chain[at_nr])] = old_state0
                        state[state_lam] = LAMIN
                        cur = numpy.array(cur0)
                        print "Start one more time", try_chain, cur, at_nr
                        break
                    #print "MO", mo 
                    chain[at_nr] = cur + mo
                    print "Rosne", at_nr, chain[at_nr], nu, i, te
                    state[tuple(chain[at_nr])] = get_site_type(i, re, la)

                    cur = chain[at_nr]
                    if i == nu-1: # the last residue of the chain, so no need to remodelling the chain 
                        succ = True
        else:
            while try_chain < 100 and not succ:
                print try_chain, succ
                for i in range(nu-1, 0, -1):
                    chain_l_dev = round(nu/125.)
                    at_nr += 1
                    if chain_l_dev * 70. <= i < chain_l_dev * 125.  and len(numpy.where(territor == int(str(te) + '1'))[0]) != 0:
                        print "ZMIENIAM", str(te) + '1'
                        open_part = numpy.where(territor == int(str(te) + '1'))
                        territor[open_part]= te
                    elif chain_l_dev * 30. <= i < chain_l_dev * 70. and len(numpy.where(territor == int(str(te) + '2'))[0]) != 0:
                        print "ZMIENIAM", str(te) + '2'
                        open_part = numpy.where(territor == int(str(te) + '2'))
                        territor[open_part]= te
                    elif chain_l_dev * 10. <= i < chain_l_dev * 30. and len(numpy.where(territor == int(str(te) + '3'))[0]) != 0:
                        print "ZMIENIAM", str(te) + '3'
                        open_part = numpy.where(territor == int(str(te) + '3'))
                        territor[open_part]= te
                    elif chain_l_dev <= i < chain_l_dev * 10. and len(numpy.where(territor == int(str(te) + '4'))[0]) != 0:
                        print "ZMIENIAM", str(te) + '4'
                        open_part = numpy.where(territor == int(str(te) + '4'))
                        territor[open_part]= te 
        
                    #print "AT_nr!!!", at_nr, i,  cur, type(cur)
                    mo, tang = rand_next(cur, state, chain, territor, nu, te, at_nr, i)
                    if tang:
                        try_chain +=1
                        assert try_chain != 100, "Chain is not possible to initialize"
                        at_nr = at_nr - nu + i
                        old_state0 = state[tuple(chain[at_nr])]
                        territor[ter_zero]= te
                        state[ter_zero] = 0.
                        territor[ter_one]= int(str(te) + '1')
                        state[ter_one] = 0.
                        territor[ter_two]= int(str(te) + '2')
                        state[ter_two] = 0.
                        territor[ter_three]= int(str(te) + '3')
                        state[ter_three] = 0.
                        territor[ter_four]= int(str(te) + '4')
                        state[ter_four] = 0.
                        state[tuple(chain[at_nr])] = old_state0
                        state[state_lam] = LAMIN
                        cur = numpy.array(cur0)
                        print "Start one more time", try_chain, cur, at_nr
                        break
                    #print "MO", mo 
                    chain[at_nr] = cur + mo
                    print at_nr, chain[at_nr], nu, i, te
                    state[tuple(chain[at_nr])] = get_site_type(i-1, re, la)

                    cur = chain[at_nr]
                    if i == 1: # the last residue of the chain, so no need to remodelling the chain 
                        succ = True
        
    mid = bound/2
    for i in range(m):
        x = random.randint(0, bound)
        y = random.randint(0, bound)
        z = random.randint(0, bound)
        tries = 0
        distance = dist_from_mi(x, y, z, mid)
        while distance > mid-3 or (not (no_collisions((x, y, z), state)) and tries < 100):
            x = random.randint(0, bound)
            y = random.randint(0, bound)
            z = random.randint(0, bound)
            tries += 1
            distance = dist_from_mi(x, y, z, mid)
        binders[i] = [x, y, z]
        state[tuple(binders[i])] = BINDER

    return chain, binders, attached_to_lamins, state


def good_neighbors(x, i, chain, s_pos_chain, l_pos_chain):
    #check neighbors (chain)
    if i > 0 and (i not in s_pos_chain):
        #print "CHECK",s_pos_chain, i
        d1 = dist(chain[i - 1], x)
        if d1 > 0 and d1 < GOOD_NEIGH:
            pass
        else:
            return False
    if i not in l_pos_chain:
        d2 = dist(chain[i + 1], x)
        if d2 > 0 and d2 < GOOD_NEIGH:
            pass
        else:
            return False
    return True

def bonds(chain, stat):
    """Return the number of bindings of the chain with binders and lamins."""
    bonds = 0
    for j in range(chain.shape[0]):
        molecule_pos = chain[j]
        molecule = stat[tuple(molecule_pos)]
        #print j
        if molecule == REGDNA:
            continue
 
        elif molecule == BSITE_R:
            binding = [BINDER]
        elif molecule == BSITE_L:
            binding = [LAMIN, BINDER]
        else: print "INNY STAN!!!", j, molecule
        #print j, molecule
        
#        one_ch_at = 0
        for bmove in BMOVES: 
            new = molecule_pos + bmove
            enc = stat[tuple(new)]
            if enc in binding:
                bonds += 1
                #one_ch_at += 1
                #if one_ch_at > 6: break # one atom can have only 6 binding partners 
    return bonds

def modify(sta_pos_chain, la_pos_chain, chain, binders, state, bound = BOUND):
    """ Try to move a binder or a bead in the chain.
    
    Return None if nothing was moved.
    Return (False, number of the moved binder, move vector) if a binder was moved.  
    Return (True, number of the moved bead in chain, move vector) if a bead in the chain was moved. 
    
    Some keywords:
    sta_pos_chain -- starting positions of chains (starting from the second).
    la_pos_chain -- last positions of chains.""" 

    # With probability 1/2 move a random binder.    
    if random.randint(0, 1): 
        i = random.randint(0, len(binders) - 1) # choose random binder
        move = random.choice(MOVES)
        new = move + binders[i]
        #print "OLD_BIN", binders[i]
        if no_collisions(tuple(new), state):
            return False, i, move
    # Else move a random residue in the chain.
    else:
        i = random.randint(0, len(chain) - 1)
        move = random.choice(MOVES)
        new = move + chain[i]
        #print "modify", i, move, chain[i], new

        # Test if there are no collisions (the same place occupied by different atoms) 
        # and no intersection of bonds. 
        if good_neighbors(new, i, chain, sta_pos_chain, la_pos_chain) and no_collisions(tuple(new), state):
            if i not in sta_pos_chain and i not in la_pos_chain: #TODO what if i == 0, not checked here
                if dist(chain[abs(i-1)], new) <= numpy.sqrt(2) and dist(chain[i+1], new) <= numpy.sqrt(2)\
                        and not intersect(new, chain[abs(i-1)], state, chain)\
                        and not intersect(new, chain[i+1], state, chain):
                    #print "Nie przecin", i
                    return True, i, move
                else: pass
                         
            elif i in la_pos_chain:
                #print "Last", i, la_pos_chain 
                if dist(chain[abs(i-1)], new) <= numpy.sqrt(2) and not intersect(new, chain[abs(i-1)], state, chain):
                #print "Nie przecin", i
                    return True, i, move
            elif i in sta_pos_chain: # TODO what if i == 0
                #print "FIRST", i, sta_pos_chain
                if dist(chain[i+1], new) <= numpy.sqrt(2) and not intersect(new, chain[i+1], state, chain):
                #print "Nie przecin", i
                    return True, i, move                
            else:
                pass
                #print "Za duza odleglosc", i
        else:
            pass 
            #print i, "No movement"
    return None    
    
DIST = 3        


def count_bonds(pos, accepted, state):
    bonds = 0
    for bmove in BMOVES:
        if state[tuple(pos + bmove)] in accepted:
            if state[tuple(pos + bmove)] == LAMIN:
                bonds += 2
            else: bonds += 1
    #       if bonds > 6: break
    #print bonds
    return bonds

def radius_gyr(chai, last_pos):
    length = chai.shape[0]
    gyr_chains = []
    for la in last_pos:
        r = 0.0
        ind = last_pos.index(la)
        if ind == 0:
            for at in range(la):
                for ato in range(at, la):
                    r = r + math.sqrt(numpy.sum((chai[at] - chai[ato])**2))
            r_gyr = 2*r/(length*(length-1))
            gyr_chains.append(r_gyr)
        else:
            for at in range(last_pos[ind-1],la):
                for ato in range(at, la):
                    r = r + math.sqrt(numpy.sum((chai[at] - chai[ato])**2))
            r_gyr = 2*r/(length*(length-1))
            gyr_chains.append(r_gyr)
    
    return gyr_chains

def write_as_pdb(ch_nr, chain, binders, attached_to_lamins, state, f, out_num, sim_step, E):
    chain_states = [state[tuple(chain[at_nr])] for at_nr in range(chain.shape[0])]
    writing_pdb.write_pdb_frame(f, ch_nr, chain, chain_states, binders, out_num, sim_step, E)

DELTA = 2
GYRATION = False
CHECK_E = False

def metropolis(nr_chrom, revers, chain, binders, attached_to_lamins, state, out_fname, n, binary,
               simulation_name, simulation_description, fpkf):
    """Run the simulation using Metropolis algorithm.

    Some keywords:
    nr_chrom - a list of lengths of chains
    revers - some lengths of chains with special meanings
    binders - coordinates of atoms of binders
    attached_to_lamins - positions of beads in the chains binding lamins
    state - 3D array of states 
    out_fname - name of file for writing
    n - number of steps
    binary - is the output binary
    num_deltas - the number of delta-frames per keyframe
    """
    
    kf_header_offsets = [] # a list of offsets of the consecutive keyframe headers
    num_deltas = fpkf - 1 # the number of delta frames is equal to the number of frames per keyframe - 1
    start_pos_chain = [] # starting positions of different subchains in the chain (starting from the second subchain)
    po = 0
    for pos in nr_chrom:
        po = po + pos
        start_pos_chain.append(po)
    
    last_pos_chain = [s-1 for s in start_pos_chain] # position of ends of different subchains
    out_file = open(out_fname, "w")
    st_nr = 0 # the number of steps made in which energy changed and pdb was outputted
    move_nr = 0 # the number of moves made (and accepted)
    kf_num = 0 # the number of the keyframe
    moves_list = []
    E = bonds(chain, state)
    print "Starting energy:", E #currently equal to the number of bonds
    if binary:
        ch_states = [state[tuple(chain[at_nr])] for at_nr in range(chain.shape[0])]
        cur_offset = writing_binary.write_bin_sim_header(out_file, simulation_name, simulation_description,
                                                         nr_chrom, ch_states, fpkf)
#        out_file, atoms_num, simulation_name, simulation_description,
#                         chr_lengths, chain, binders, state)
#        cur_offset = write_as_binary(nr_chrom, chain, binders, out_file, st_nr, 0, cur_offset, name + ";bonds=" + str(E))
    else:
        write_as_pdb(nr_chrom, chain, binders, attached_to_lamins, state, out_file, st_nr, 0, E)
    for step in range(n):
        resp = modify(start_pos_chain, last_pos_chain, chain, binders, state)
        ch = numpy.array(chain, copy = True)
        b = numpy.array(binders, copy = True)
        if CHECK_E:
            st = numpy.array(state,copy=True)
        else: pass #TODO unnecessary

        if resp: # If a binder or a bead in the chain was moved.
            in_chain, i, move = resp
            if in_chain: # If a bead in the chain was moved.
                #print i, "change"
                old = numpy.copy(ch[i])
                ch[i] = ch[i] + move
                if CHECK_E:
                    st[tuple(ch[i])] = st[tuple(old)]
                    st[tuple(old)] = EMPTY
                else: pass
                if state[tuple(old)] == BSITE_R:
                    Enew = E + count_bonds(ch[i], [BINDER], state) - count_bonds(old, [BINDER], state)
                    if CHECK_E:
                        Eslow = bonds(ch, st)
                        if Enew != Eslow:
                            print "R", 'Enew', Enew, 'Eslow', Eslow
                    else: pass
                elif state[tuple(old)] == BSITE_L:
                    Enew = E + count_bonds(ch[i], [LAMIN, BINDER], state) - count_bonds(old, [LAMIN, BINDER], state)
                    if CHECK_E:
                        Eslow = bonds(ch, st)
                        if Enew != Eslow:
                            print "L", 'Enew', Enew, 'Eslow', Eslow
                    else: pass
                    if tuple(ch[i]) not in attached_to_lamins and count_bonds(ch[i], [LAMIN], state) > 0:
                        #print "NOT i > 0", tuple(ch[i]),  attached_to_lamins, count_bonds(ch[i], [LAMIN], state), ch[i]
                        attached_to_lamins.append(tuple(ch[i]))
                    elif tuple(ch[i]) in attached_to_lamins and count_bonds(ch[i], [LAMIN], state) == 0:
                        #print "IN i ==0", tuple(ch[i]),  attached_to_lamins, count_bonds(ch[i], [LAMIN], state), ch[i]
                        attached_to_lamins.remove(ch[i])
                    #elif tuple(ch[i]) in attached_to_lamins and count_bonds(ch[i], [LAMIN], state) > 0:
                    #    print "JEST, ale ma lamine!", ch[i], count_bonds(ch[i], [LAMIN], state)
                    else:
                        pass
                        #print "WYJEATEK ", tuple(ch[i]), attached_to_lamins, count_bonds(ch[i], [LAMIN], state)
                else: # REGDNA
                    Enew = E
            else: # If a binder was moved.
                old = numpy.copy(b[i])
                b[i] = b[i] + move
                if CHECK_E:
                    st[tuple(b[i])] = st[tuple(old)]
                    st[tuple(old)] = EMPTY
                else: pass
                Enew = E + count_bonds(b[i], [BSITE_R, BSITE_L], state) - count_bonds(old, [BSITE_R, BSITE_L], state)
                if CHECK_E:
                    Eslow = bonds(ch, st)
                    if Enew != Eslow:
                        print "B", 'Enew', Enew, 'Eslow', Eslow
                else: pass
        else: # Nothing was moved.
            Enew = E

        if Enew >= E or random.uniform(0.0, 1.0) < math.exp((Enew - E) * DELTA): #accept the step
            chain = ch
            binders = b
            if resp: # If something was moved.
                state[tuple(old + move)] = state[tuple(old)]
                state[tuple(old)] = EMPTY
                if binary:
                    move_nr += 1
                    moves_list.append((old, move, step))
                    if move_nr % num_deltas == 0:
                        assert len(moves_list) % num_deltas == 0
                        kf_header_offsets.append(cur_offset)
                        cur_offset = writing_binary.write_as_binary(nr_chrom, chain,
                            binders, out_file, moves_list, step, kf_num, cur_offset, Enew, 1)
                        kf_num += 1
#                    (ch_lengths, chain, binders, f, move_nr, moves_list, metr_step, kf_num, cur_offset, name="chromosome and binders")
                        moves_list = []
            if E != Enew: # if the energy changed
                E = Enew
                st_nr += 1
                if GYRATION:
                    print "iter", step, "step", st_nr, "energy:", E, "R_gyr ", radius_gyr(chain, last_pos_chain)
                else:
                    print "iter", step, "step", st_nr, "energy:", E
                if not binary:
                    write_as_pdb(nr_chrom, chain, binders, attached_to_lamins, state, out_file,
                                 st_nr, step, E)

    if binary:
        writing_binary.set_kf_header_offsets(out_file, kf_header_offsets)
        writing_binary.complete_bin_sim_header(out_file, kf_num, fpkf, cur_offset)

    # dump the last state to the pickle
    l_obj = [chain, binders, attached_to_lamins, state, nr_chrom, revers]
    pickle_fname = out_fname.split('.')[0] + ".pick"
    pickle_file = open(pickle_fname, 'w')
    pickle.dump(l_obj, pickle_file)
    out_file.close()

def output_name(ou, m, n, binary):
    if ou == '':
        if type(m) == list:
            f_n = "MC_traj_%ibin_%ichain" % (len(m), len(n))
        else:
            f_n = "MC_traj_%ibin_%ichain" % (m, n)
    elif "." in ou:
        f_n = ou.split('.')[0]
    else: 
        f_n = ou
    if binary:
        f_n = ou + ".bin"
    else:
        f_n = ou + ".pdb"
    return f_n

args = pars_inp()

if args.In_str == '':
    rand_init = True
else: 
    rand_init = False

t1 = time.time()
if rand_init: 
    N = args.Ch_lenght.split(',') # length of the chains
    N = map(int, N)
    R = args.Revers.split(',')
    R = map(int, R)
    M = args.Binders # nr of binders
    fn = output_name(args.Out_str, M, N, args.binary)
    c, b, a, state = initialize_random(N, M, fn, R, args.regularBSitesFile, args.laminBSitesFile, args.binary)
else:
    c, b, a, state, N, R = initialize_import(args.In_str)
    M = b.shape[0]
    fn = output_name(args.Out_str, M, N, args.binary)

print "The lenght of the chain is ", N, ", the number of binders is ", M, ", the nucleus radius is ", R, ", the number of steps is ", args.Steps, "random_seed ", #ran_seed
t2 = time.time()
print "initialization: ", t2 - t1
BOUND = numpy.max(c)

a = []
metropolis(N, R, c, b, a, state, fn, args.Steps, args.binary, args.simulation_name, args.simulation_description, args.fpkf)
t1 = t2
t2 = time.time()
print "metropolis: ", t2 - t1
