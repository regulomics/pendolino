"""A module with functions for converting from the pdb to the binary format."""
import constants
from constants import DIST
import writing_binary
import numpy
from collections import namedtuple

def iter_chains(pdb_f):
    """Yield the next line from a pdb_f file if it describes an atom in a chain."""
    line = pdb_f.readline()
    while len(line) >= 4 and line[0:4] == 'ATOM' and line[13] == 'C':
        yield line
        line = pdb_f.readline()

def iter_chain(pdb_f, line, terminal):
    """Yield a line with an atom."""
    while line[:4] == 'ATOM':
        yield line
        line = pdb_f.readline()
    assert line[0:3] == terminal

def read_chain_states_from_pdb(pdb_f):
    """Return a list of chain lengths, states, and the number of binders."""
    line1 = pdb_f.readline()
    assert line1[:6] == "HEADER"
    line2 = pdb_f.readline()
    assert line2[:5] == "TITLE"
    ch_lens = []
    ch_states = []
    names_regions = {name: reg for reg, name in constants.pdb_reg_names.items()}
    for line in iter_chains(pdb_f):
        ch_len = 0
        for line2 in iter_chain(pdb_f, line, 'TER'):
            ch_len += 1
            ch_states.append(names_regions[line2[17:20]])
        ch_lens.append(ch_len)    
    num_binders = 0
    for _ in iter_chain(pdb_f, line2, 'CON'):
        num_binders += 1
    return ch_lens, ch_states, num_binders

def line_to_coord(line):
    """Return coordinates from the pdb line."""
    return numpy.array([int(float(line[-36 + 8*i:-28 + 8*i].lstrip())/DIST) for i in range(3)])

def parse_pdb_frame(pdb_f, ch_len, num_binders, line):
    """Parse pdb frame. Return step, energy, and arrays of coordinates of atoms in chains and coordinates of binders."""
    step_str = "step "
    step = int(line[line.find(step_str) + len(step_str):])
    line = pdb_f.readline()
    bds = "bonds="
    E = int(line[line.find(bds) + len(bds):])
    at_nr = 0
    chain = numpy.zeros((ch_len, 3), dtype = numpy.int)

    for line in iter_chains(pdb_f):
        for line2 in iter_chain(pdb_f, line, 'TER'):
            chain[at_nr] = line_to_coord(line2)
            at_nr += 1

    binders = numpy.zeros((num_binders, 3), dtype = numpy.int)
    at_nr = 0

    for line in iter_chain(pdb_f, line2, 'CON'):
            binders[at_nr] = line_to_coord(line)
            at_nr += 1
    return step, E, chain, binders

def get_new_pdb_frame_line(in_file):
    """Return the line of the beginning of a new pdb frame if such a frame exists and an empty string otherwise."""
    line = in_file.readline() 
    while line:
        if len(line) >= 6 and line[0:6] == "HEADER":
            return line
        line = in_file.readline()     
    return ""                     

def get_pos_changes_info(chain, chainNew):
    """Return a list of pairs (i, v) where v is the displacement and i is its index."""
    chain_changes = []
    for i, pos in enumerate(chain):
        if any(pos != chainNew[i]):
            chain_changes.append((i, chainNew[i] - pos))
    return chain_changes

def iter_moves(change):
    """Decompose change into small changes moving at most two coordinates by one. Yield the small changes."""
    pairs = [(0, 1), (0, 2), (1, 2)]
    for i, j in pairs:
        while (change[i] != 0) and (change[j] != 0):
            move = numpy.zeros(3, dtype = numpy.int)
            move[i] = numpy.sign(change[i])
            move[j] = numpy.sign(change[j])
            change -= move
            yield move
    for i in range(3):
        while change[i] != 0:
            move = numpy.zeros(3, dtype = numpy.int)  
            move[i] = numpy.sign(change[i])
            change -= move
            yield move
                        
def iter_moves_chain(moves_list, chain, chain_changes, step, num_deltas):
    """Yield a list of up to num_deltas triples: (previous position of an atom, displacement, step) in chain."""
    for i, change in chain_changes:
        for move in iter_moves(change):
            moves_list.append((chain[i], move, step))
            chain[i] += move
            if len(moves_list) == num_deltas:
                yield moves_list
                moves_list = []
    yield  moves_list

def iter_data_for_bin(chain, binders, chain_changes, binders_changes, step, num_deltas):
    """Yield a list of up to num_deltas triples: (previous position of atom, displacement, step) in chain or binders.
    
    Also move the chains and binders by the returned moves."""
    moves_list = []
    for moves_list in iter_moves_chain(moves_list, chain, chain_changes, step, num_deltas):
        if len(moves_list) == num_deltas:
            yield moves_list
    for moves_list in iter_moves_chain(moves_list, binders, binders_changes, step, num_deltas):
        if len(moves_list) == num_deltas:
            yield moves_list
    if moves_list: # there are still some moves left, fill moves_list with zeros
        z = numpy.zeros(3, dtype = numpy.int)
        for _ in range(num_deltas - len(moves_list)):
            moves_list.append((z, z, step))
        yield moves_list

def convert_lamin_pdb_to_bin(in_lamin, out_lamin):
    binders_list = []
    with open(in_lamin, 'r') as in_file:
        line = in_file.readline()
        assert line == "HEADER LAMINA\n"
        line = in_file.readline()
        while line:
            binders_list.append(list(line_to_coord(line)))            
            line = in_file.readline()
        writing_binary.write_lamins_binary(out_lamin, binders_list)

def convert_sim_pdb_to_bin(input_name, out_name, fpkf):
    num_deltas = fpkf - 1
    with open(input_name, 'r') as in_file, open(out_name, 'w') as out_file:
        ch_lengths, ch_states, num_binders = read_chain_states_from_pdb(in_file)
        cur_offset = writing_binary.write_bin_sim_header(out_file, "", "",
                                                         ch_lengths, ch_states, fpkf)
        in_file.seek(0)
        kf_num = 0
        kf_header_offsets = [] # A list of offsets of the consecutive keyframe headers. 
        line =  in_file.readline()
        step, E, chain, binders = parse_pdb_frame(in_file, sum(ch_lengths), num_binders, line)
        line = get_new_pdb_frame_line(in_file)

        while line: # while there are still some frames and line is the first line of o 
            step, E, new_chain, new_binders = parse_pdb_frame(in_file, sum(ch_lengths), num_binders, line)
            chain_changes = get_pos_changes_info(chain, new_chain)
            binders_changes = get_pos_changes_info(binders, new_binders)
            for moves_list in iter_data_for_bin(chain, binders, chain_changes, 
                                                binders_changes, step, num_deltas):
                kf_header_offsets.append(cur_offset)
                cur_offset = writing_binary.write_as_binary(ch_lengths, chain, binders, out_file,
                    moves_list, step, kf_num, cur_offset, E, 1)
                kf_num += 1
            line = get_new_pdb_frame_line(in_file)
        writing_binary.set_kf_header_offsets(out_file, kf_header_offsets)
        writing_binary.complete_bin_sim_header(out_file, kf_num, fpkf, cur_offset)


def convert_pdb_to_bin(input_name, out_name, fpkf, in_lamin, out_lamin):
    """Convert a pdb files to a binary ones. 
    
    Keywords:
    input_name - name of the input pdb simulation file
    out_name - name of the output binary simulation file
    fpkf - the number of frames per keyframe"""

    if in_lamin:
        convert_lamin_pdb_to_bin(in_lamin, out_lamin)
    convert_sim_pdb_to_bin(input_name, out_name, fpkf)
    
    
    
