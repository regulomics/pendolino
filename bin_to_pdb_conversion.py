"""A module with functions for converting from the binary to the pdb format."""

import constants, numpy, struct, message_pb2, writing_pdb
from collections import namedtuple

def parse_h_data(h_data):
    """Parse header data. Return a list of lengths of chains and a list of states of beads in chains."""
    ch_lens = []
    ch_states = []
    assert h_data.binders_types_count == 2
    evec_state = {v: i for i, v in constants.binding_forces.iteritems()} # a map energy vector -> state

    for chain in h_data.chains:
        ch_len = 0
        for bead in chain.beads:
            ch_len += 1
            e_v = [0 for _ in range(h_data.binders_types_count)]
            for vector in bead.energy_vector:
                e_v[vector.binder_type] = vector.force
            ch_states.append(evec_state[(tuple(e_v))])
        ch_lens.append(ch_len)
    return ch_lens, ch_states

def get_binders_list(kf_data):
    binders_list = []
    for binder in kf_data.binders:
        p = binder.position
        binders_list.append([p.x, p.y, p.z])
    return binders_list


def write_kf_as_pdb(pdb_f, kf_data, ch_lens, ch_states, out_num, E):
    """Write keyframe as pdb.
    
    Some keywords:
    pdb_f - a pdb file
    kf_data - keyframe data
    ch_lens - a list of lengths of chains
    ch_states - a chain of states of beads in chains
    out_num - the number of the converted keyframe (i.e. among ones in which energy changed)"""

    chain = numpy.zeros((len(ch_states), 3), dtype = numpy.int)
    at_nr = 0
    for ch in kf_data.chains:
        for bead in ch.bead_positions:
            chain[at_nr] = numpy.array([bead.x, bead.y, bead.z])
            at_nr += 1
    binders = numpy.array(get_binders_list(kf_data))
    writing_pdb.write_pdb_frame(pdb_f, ch_lens, chain, ch_states, binders, out_num, kf_data.step_counter, E)

class Kfh_parser:
    """A class used for parsing keyframe headers."""
    def __init__(self):
        list_height = 10
        self.kfh_fs = ">" + (2 + list_height)*"Q" + "L"
        self.kfh_len = (2 + list_height)*8 + 4
        kf_fields = "kf_num delta_list_offset "
        for i in range(list_height):
            kf_fields += "e" + str(i) + " "
        kf_fields += "kf_data_size"
        self.kfh_nt = namedtuple("kf_header", kf_fields)
        self.max_step = 2**(list_height - 1)
        
    def parse(self, bin_f):
        kfh_str = bin_f.read(self.kfh_len)
        t = struct.unpack_from(self.kfh_fs, kfh_str)
        return self.kfh_nt._make(t)

def move_frame(bin_f, step, kfh, kfh_p):
    """Go step frames ahead."""
    assert step > 0
    end = False
    while not end:
        k = 0
        i = 1
        while 2*i <= min(step, kfh_p.max_step):
            i = i*2
            k = k + 1
        bin_f.seek(kfh.__getattribute__("e" + str(k)))
        step -= i
        if step == 0:
            end = True
        else:
            kfh = kfh_p.parse(bin_f)

def parse_bin_header(bin_f):
        """Parse the header of the binary file."""
        header_str = bin_f.read(56)
        fs = ">8sQQQQQL"
        ph_nt = namedtuple("proto_header","protostr file_size hdr_pos keyframe_pos keyframe_no frame_no frames_per_key")
        t = struct.unpack_from(fs, header_str)
        ph = ph_nt._make(t)
        h_data_len = ph.keyframe_pos - ph.hdr_pos
        h_data_str = bin_f.read(h_data_len)
        h_data = message_pb2.Header()
        h_data.ParseFromString(h_data_str)
        ch_lens, ch_states = parse_h_data(h_data)
        return ph.keyframe_no, ch_lens, ch_states

def get_kf_data(bin_f, kfh):
        kf_data = message_pb2.Keyframe()
        kf_data_str = bin_f.read(kfh.kf_data_size)
        kf_data.ParseFromString(kf_data_str)
        return kf_data

def convert_lamin_bin_to_pdb(in_lamin, out_lamin):
    with open(in_lamin, 'r') as bin_f:
        parse_bin_header(bin_f)
        kfh_p = Kfh_parser()
        kfh = kfh_p.parse(bin_f)
        kf_data = get_kf_data(bin_f, kfh)
        writing_pdb.write_lamins_pdb(out_lamin, get_binders_list(kf_data))
        
def convert_sim_bin_to_pdb(input_name, out_name, step):
    with open(input_name, 'r') as bin_f, open(out_name, 'w') as pdb_f:
        keyframe_no, ch_lens, ch_states = parse_bin_header(bin_f)
        E = -1
        out_num = 0 
        i = 0
        kfh_p = Kfh_parser()
        assert step >= 0
        while True:
            kfh = kfh_p.parse(bin_f)
            kf_data =  get_kf_data(bin_f, kfh)
            E_new = kf_data.callbacks[0].int_value
            if (step > 0) or (step == 0 and E_new != E):
                E = E_new
                write_kf_as_pdb(pdb_f, kf_data, ch_lens, ch_states, out_num, E)
                out_num += 1
            if step == 0:
                i += 1
            else:                     
                i += step
            if i >= keyframe_no:
                break
            if step == 0 or step == 1:
                bin_f.seek(kfh.e0)
            else:
                move_frame(bin_f, step, kfh, kfh_p)

def convert_bin_to_pdb(input_name, out_name, step, in_lamin, out_lamin):
    if in_lamin:
        convert_lamin_bin_to_pdb(in_lamin, out_lamin)
    convert_sim_bin_to_pdb(input_name, out_name, step)
