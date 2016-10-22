"""A module with functions for outputting to binary format."""

import message_pb2, constants, struct, numpy
def get_h_data_str(simulation_name, simulation_description,
                         chr_lengths, ch_states):
    """Produce a string containing the header data."""

    h_data = message_pb2.Header()
    h_data.simulation_name = simulation_name
    h_data.simulation_description = simulation_description
    h_data.binders_types_count = 2
    h_data.binder_types_names.extend(["Lamin", "RegularBinder"])

    atom_nr = 0
    for chr_n, chr_l in enumerate(chr_lengths):
        proto_chain = h_data.chains.add()
        proto_chain.chain_name = constants.chain_list[chr_n]
        for _  in range(chr_l):
            bead = proto_chain.beads.add()
            for i, b_force in enumerate(constants.binding_forces[ch_states[atom_nr]]):
                bead.energy_vector.add(binder_type=i, force=b_force)
            atom_nr += 1
    return h_data.SerializeToString()

def write_bin_sim_header(out_file, simulation_name, simulation_description,
                         ch_lens, ch_states, fpkf):
    """Return the offset of the first keyframe header.

    Write the protostream header and header data of the binary file with simulation data.
    Keywords:
    out_file - the output binary file
    simulation_name - the name of the simulation
    simulation_description - the description of the simulation
    ch_lengths - a list of lengths of all chains
    ch_states - a list of states of beads in chains
    fpkf - the number of frames per keyframe"""

    h_data_str = get_h_data_str(simulation_name, simulation_description,
                                ch_lens, ch_states)
    fs = ">8sQQQQQL"
    h_data_pos = 56
    keyframe_h_pos = h_data_pos + len(h_data_str)
    header = struct.pack(fs, "PROTOSTR", 0, h_data_pos, keyframe_h_pos, 0, 0, fpkf)
    out_file.write(header)
    out_file.seek(56)
    out_file.write(h_data_str)
    return keyframe_h_pos

def get_kf_data_str(ch_lengths, chain, binders, binder_type, metr_step, Energy):
    """Return a keyframe data string"""
    kf_data = message_pb2.Keyframe()
    kf_data.step_counter = metr_step
    for i in range(binders.shape[0]):
        pos = binders[i]
        binder = kf_data.binders.add()
        binder.binder_type = binder_type
        binder.position.x = pos[0]
        binder.position.y = pos[1]
        binder.position.z = pos[2]
    atom_nr = 0
    kf_data.callbacks.add(int_value = Energy) 
    for ch_len in ch_lengths:
        proto_chain = kf_data.chains.add()
        for _ in range(ch_len):
            pos = chain[atom_nr]
            proto_chain.bead_positions.add(x=pos[0], y=pos[1], z=pos[2])
            atom_nr += 1
    return kf_data.SerializeToString()

def get_deltas_str(moves_list):
    """Return a string with delta-frames."""
    deltas_str = ""
    for old, move, step in moves_list:
        delta = message_pb2.Delta()
        delta.__getattribute__("from").x = old[0]
        delta.__getattribute__("from").y = old[1]
        delta.__getattribute__("from").z = old[2]
        delta.disp.x = move[0]
        delta.disp.y = move[1]
        delta.disp.z = move[2]
        delta.step_counter = step
        deltas_str += delta.SerializeToString() 
    return deltas_str

def get_kf_header_str(cur_offset, kf_data_size, kf_num):
    """Return a pair: a string with a keyframe header and a delta list offset."""   
    list_height = 10
    fs = ">" + (2 + list_height)*"Q" + "L"
    k_header_len = (2 + list_height)*8 + 4
    delta_list_offset = cur_offset + k_header_len + kf_data_size #TODO
    kf_header_str = struct.pack(fs, kf_num, delta_list_offset, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, kf_data_size)
    return kf_header_str, delta_list_offset

def write_as_binary(ch_lengths, chain, binders, f, moves_list, metr_step, kf_num, cur_offset, E, binder_type):
    """Write a keyframe (i.e. keyframe header, keyframe_data, and delta-frames) and return the current offset of the end of the binary file.

    Some keywords:
    ch_lengths - a list of lengths of chains
    chain - a list of coordinates of atoms in chains
    binders - a list of coordinates of binders
    move_nr - the number of the move
    metr_step - the number of steps made
    cur_offset - current offset, equal to the current end of file f
    moves_list - a list of triples (old, move, step), where old is the old position of an atom,
    move is its displacement, and step is the number of the step of the movement
    E - the energy"""

    kf_data_str = get_kf_data_str(ch_lengths, chain, binders, binder_type, metr_step, E)
    kf_data_size = len(kf_data_str)    
    kf_header_str, delta_list_offset = get_kf_header_str(cur_offset, kf_data_size, kf_num)
    f.write(kf_header_str)
    f.write(kf_data_str)
    deltas_str = get_deltas_str(moves_list)
    f.write(deltas_str)
    cur_offset = delta_list_offset + len(deltas_str)
    return cur_offset

def set_kf_header_offsets(out_file, kf_header_offsets):
    """Set offsets of keyframes in keyframe headers.
    
    Some keywords:
    kf_header_offsets - a list of offsets of all keyframe headers in the binary file."""
    list_height = 10
    num_kf = len(kf_header_offsets)
    fs = ">Q"
    for i, offset in enumerate(kf_header_offsets):
        k = 0
        j = i + 1
        out_file.seek(offset + 16) 
        while j < num_kf and k < list_height:
            kf_pos_str = struct.pack(fs, kf_header_offsets[j])
            out_file.write(kf_pos_str)
            k = k + 1
            j = i + 2**k

def complete_bin_sim_header(f, kf_num, fpkf, cur_offset):
    """Complete the information in the header of the binary file."""
    fs = ">Q"
    f.seek(8)
    file_size_str = struct.pack(fs, cur_offset)
    f.write(file_size_str)
    f.seek(32)
    kf_num_str = struct.pack(fs, kf_num)
    f.write(kf_num_str)
    f.seek(40)
    frames_num_str = struct.pack(fs, kf_num*fpkf)
    f.write(frames_num_str)

def write_lamins_binary(lam_name, binders_list):
    with open(lam_name, "w") as save_lam:
        cur_offset = write_bin_sim_header(save_lam, "", "", [], [], 1)
        cur_offset = write_as_binary([], [], numpy.array(binders_list), save_lam, [], 0, 0, cur_offset, 0, 0)
        complete_bin_sim_header(save_lam, 1, 1, cur_offset)
