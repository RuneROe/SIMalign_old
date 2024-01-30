import numpy as np
from pymol import cmd
from scipy.spatial import cKDTree
from Bio import AlignIO


def resi_to_index(resi,ref_structure,structure_list,align):
    model = cmd.get_model(ref_structure+" and name CA and chain A")
    resi = resi - model.atom[0].resi + 1
    ref_index = structure_list.index(ref_structure)
    count = 0
    for index, AA in enumerate(align[ref_index].seq):
        if AA != "-":
            count += 1
        if count == int(resi):
            return index

def index_to_resi(index,ref_structure,structure_list,align):
    ref_index = structure_list.index(ref_structure)
    model = cmd.get_model(ref_structure+" and name CA and chain A")
    resi = model.atom[0].resi - 1
    for i, AA in enumerate(align[ref_index].seq):
        if AA != "-":
            resi += 1
        if i == index:
            return resi

def bigger_AA(target_AA,ref_AA):
    if ref_AA == target_AA:
        return False
    elif ref_AA == "G":
        return False
    elif ref_AA == "A" and target_AA not in {"G","P"}:
        return False
    elif ref_AA == "V" and target_AA == "I":
        return False
    elif ref_AA == "F" and target_AA == "Y":
        return False
    elif ref_AA == "L" and target_AA in {"F","Y"}:
        return False
    elif ref_AA == "S" and target_AA in {"C","T"}:
        return False
    else:
        return True

def dist_points(coord1, coord2):
    coord1 = np.array(coord1)
    coord2 = np.array(coord2)
    distance = np.linalg.norm(coord2 - coord1)
    return distance

def add_hotspot(closeAA_list,atom,i,structure_list,align,models):
    amino_acid_translation = {
    'ALA': 'A',
    'ARG': 'R',
    'ASN': 'N',
    'ASP': 'D',
    'CYS': 'C',
    'GLU': 'E',
    'GLN': 'Q',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'THR': 'T',
    'TRP': 'W',
    'TYR': 'Y',    
    'VAL': 'V'}
    index = None
    resi_list = []
    ref_index = resi_to_index(int(atom.resi),structure_list[i],structure_list,align)
    for j, seq in enumerate(align):
        align_char = seq[ref_index]
        try:
            if align_char != amino_acid_translation[atom.resn] and align_char != "-" and dist_points(atom.coord,models[j].atom[index_to_resi(ref_index,structure_list[j],structure_list,align)-1].coord) < 1:
                flag = True
                if not isinstance(closeAA_list[j][ref_index], str):
                    for closeAA in closeAA_list[j][ref_index]:
                        close_index = resi_to_index(closeAA,structure_list[j],structure_list,align)
                        if close_index != None:
                            if bigger_AA(seq[close_index],align[i][close_index]):
                                flag = False     
                if flag:
                    index = int(atom.resi)
                    resi_list.append(index_to_resi(ref_index,structure_list[j],structure_list,align))
                else:
                    resi_list.append("-")
            else:
                resi_list.append("-")
        except:
            resi_list.append("-")
    return index, resi_list


def discard_surface_residues(hotspot,structure):
    import findsurfaceatoms
    exposed_residues = findsurfaceatoms.findSurfaceResidues(selection=structure+" and chain A and not HETATM")
    exposed_set = set()
    for resi in exposed_residues:
        exposed_set.add(resi[1])
    new_hotspot = {}
    for k,v in hotspot.items():
        if k not in exposed_set:
            new_hotspot[k] = v
    return new_hotspot, exposed_set

def get_close_aa(close_AAs,modelatoms,kd,atom,resi):
    tmp_set = set([int(modelatoms[x].resi) for x in kd.query_ball_point(atom.coord,4)])
    tmp_set.discard(resi)
    return close_AAs.union(tmp_set)

def get_close_aa_list(structure_list,align):
    closeAA_list = []
    for i, seq in enumerate(align):
        closeAA_list.append([])
        for ele in seq:
            closeAA_list[i].append(ele)

    for i, structure in enumerate(structure_list):
        model = cmd.get_model(structure+" and (not name CA and not name N and not name C and not name O) and chain A and not HETATM")
        kd = cKDTree([atom.coord for atom in model.atom])
        resi = 0
        close_AAs = set()
        for atom in model.atom:
            if resi != int(atom.resi):
                if resi != 0:
                    # if resi_to_index(resi,structure,structure_list,align) != None:
                    closeAA_list[i][resi_to_index(resi,structure,structure_list,align)] = close_AAs
                    # else:
                        # pass
                    close_AAs = set()
                resi = int(atom.resi)
                close_AAs = get_close_aa(close_AAs,model.atom,kd,atom,resi)
            else:
                close_AAs = get_close_aa(close_AAs,model.atom,kd,atom,resi)
        # if resi_to_index(resi,structure,structure_list,align) != None:
        closeAA_list[i][resi_to_index(resi,structure,structure_list,align)] = close_AAs
        # else:
            # pass
    return closeAA_list


def run(structure_list,alignment_file_name,discard_exposed):
    print("Finding hotspots...")
    align = AlignIO.read(alignment_file_name,"clustal")


    models = []
    for structure in structure_list:
        models.append(cmd.get_model(structure+" and name CA and chain A and not HETATM"))

    closeAA_list = get_close_aa_list(structure_list,align)
    hotspot_list = []
    exposed_list = []
    for i, structure in enumerate(structure_list):
        hotspot = {}
        for atom in cmd.get_model(structure+" and name CA and not HETATM and chain A").atom:
            k,v = add_hotspot(closeAA_list,atom,i,structure_list,align,models)
            if k != None:
                hotspot[k] = v

        if discard_exposed:
            hotspot, exposed_set = discard_surface_residues(hotspot,structure)
            exposed_list.append(exposed_set)
        else:
            exposed_set = None
        hotspot_list.append(hotspot)
    return hotspot_list, exposed_list


def print_hotspot(hotspot_list,structure_list,print_hospots_from_structure):
    model_list = []
    for structure in structure_list:
        model_list.append(cmd.get_model(structure+" and name CA and chain A and not HETATM").atom)
    for i, hotspot in enumerate(hotspot_list):
        print("\tPrinting possible single mutations in "+structure_list[i])
        for k,v in hotspot.items():
            for j, resi in enumerate(v):
                if resi != "-":
                    print("\t\t"+model_list[i][k-1].resn+model_list[i][k-1].resi+" -> "+model_list[j][resi-1].resn+" as structure: "+structure_list[j])
        if print_hospots_from_structure == "ref_structure":
            break