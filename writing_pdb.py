"""A module with functions for outputting to the pdb format."""
import constants
from constants import DIST

def write_lamins_pdb(lam_name, binders_list):
    with open(lam_name, "w") as save_lam:
        save_lam.write("HEADER LAMINA")
        at_nr = 1
        for x, y, z in binders_list:
            line = "\nATOM  " + str(at_nr).rjust(5) + " " + "P".center(4) + " " + \
                "LAM" + "  " + str(at_nr).rjust(4) + "    " + \
                str(round(x * DIST, 3)).rjust(8) + str(round(y * DIST, 3)).rjust(8) + \
                str(round(z * DIST, 3)).rjust(8) + "  0.00 00.00"
            save_lam.write(line)
            at_nr += 1

def write_pdb_frame(f, ch_lens, chain, chain_states, binders, out_num, sim_step, E):
    """Write a single pdb frame.

    Keyword arguments:
    f - the file used for output
    ch_lens - an array of lengths of chains
    chain - a numpy array of atoms' coordinates
    chain_states - an array of states of atoms in chains
    binders - a numpy array of binders' coordinates
    out_num - the number of the outputted frame 
    sim_step  - the step number of the simulation
    E - the energy"""

    l = chain.shape[0] # the length of chain (sum of lengths of all subchains)
    n = binders.shape[0] # the number of binders
    f.write("HEADER %i %d step %i\nTITLE %s" % (out_num, l + n, sim_step, "chromosome;bonds=" + str(E)))
    at_nr = 0

    def pdb_line(at_name, at_nr, res_nr, desc, pos, chain_n):
        return "\nATOM  " + str(at_nr).rjust(5) + " " + at_name.center(4) + " " + desc \
            + " " + chain_n + str(res_nr).rjust(4) + "    " + str(round(pos[0] * DIST, 3)).rjust(8)\
            + str(round(pos[1] * DIST, 3)).rjust(8) + str(round(pos[2] * DIST, 3)).rjust(8) + "  0.00 00.00"

    # Write the positions of beads in the chain.          
    for l, ch_n in zip(ch_lens, constants.chain_list):
#        if l not in rev:
        res_nu = 0
#       else: res_nu = l + 1
        for i in range(l):
            at_n = 'C'
            r = constants.pdb_reg_names[chain_states[at_nr]]
            at_nr += 1
#            if l not in rev:
            res_nu += 1
#            else: res_nu -= 1
            f.write(pdb_line(at_n, at_nr, res_nu, r, chain[at_nr-1], ch_n))
        f.write("\nTER")            

    # Write the positions of binders. 
    chain_at = at_nr
    res_nu = 0
    for i in range(n):
        res_nu += 1
        at_nr += 1
        r = "BIN"
        at_n = 'O'
        f.write(pdb_line(at_n, at_nr, res_nu, r, binders[i], "0"))

    # Write the connections. 
    ind = 0
    sum = ch_lens[ind]
    for i in range(1, chain_at):
        #print i, sum
        if i == sum:
            ind +=1
            sum = sum + ch_lens[ind]
            continue
        else:
            line = "\nCONECT" + str(i).rjust(5) + str(i + 1).rjust(5)
        f.write(line) 
    f.write("\nEND   \n")
    
