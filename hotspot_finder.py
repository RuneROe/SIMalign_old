import numpy as np
from pymol import cmd
from scipy.spatial import cKDTree
from Bio import AlignIO


def resi_to_index(residue,align_seq,atomsCA):
    # atomsCA = cmd.get_model(ref_structure+" and name CA").atom
    # ref_index = structure_list.index(ref_structure)
    # print("ref_index",ref_index)
    count = 0
    for index, AA in enumerate(align_seq.seq):
        if AA != "-":
            if int(atomsCA[count].resi) == residue:
                return index
            count += 1

def index_to_resi(index,align_seq,atomsCA):
    # atomsCA = cmd.get_model(ref_structure+" and name CA").atom
    # ref_index = structure_list.index(ref_structure)
    count = 0
    for i, AA in enumerate(align_seq.seq):
        if AA != "-":
            if i == index:
                return int(atomsCA[count].resi)
            count += 1
        elif i == index:
            return "-"



# def resi_to_index(residue,ref_structure,structure_list,align):
#     atomsCA = cmd.get_model(ref_structure+" and name CA").atom
#     # print("r2i, resi",resi)
#     # print("r2i, int",int(model.atom[0].resi))
#     # residue = residue - int(model.atom[0].resi) + 1
#     ref_index = structure_list.index(ref_structure)
#     count = 0
#     for index, AA in enumerate(align[ref_index].seq):
#         if AA != "-":
            
#             if int(atomsCA[count].resi) == residue:
#                 return index
#         #     count += 1
#             count += 1
#         # if count == int(residue):
            

# def index_to_resi(index,ref_structure,structure_list,align):
#     atomsCA = cmd.get_model(ref_structure+" and name CA").atom
#     ref_index = structure_list.index(ref_structure)
#     # residue = int(atomsCA[0].resi) - 1
#     count = 0
#     for i, AA in enumerate(align[ref_index].seq):
#         # if AA != "-":
#         #     resi += 1
#         if AA != "-":
#             if i == index:
#                 return int(atomsCA[count].resi)
#             count += 1
#         elif i == index:
#             return "-"
#         # if i == index:
#         #     if AA == "-":
#         #         return "-"
#         #     else:
#         #         return resi

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
 
def add_hotspot(closeAA_list,structure_list,align,atoms_list,resn,resi,coord,count):
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
    residue = None
    resi_list = []
    ref_index = resi_to_index(resi,align[0],atoms_list[0])   
    # align_char = align[0][ref_index]
    # print(ref_index,resi,structure_list[0],structure_list,align,atoms_list[0])
    for j, seq in enumerate(align[1:]):
        k = j+1
    # seq = align[0]
        align_char = seq[ref_index]
        # try:
        if align_char != "-" and align_char != amino_acid_translation[resn]:
            for x in atoms_list[k]:
                if int(x.resi) == index_to_resi(ref_index,seq,atoms_list[k]):
                    break
            if dist_points(coord,x.coord) < 1:
                flag = True
                # if not isinstance(closeAA_list[j][ref_index], str):  # Check for glycine
                for closeAA in closeAA_list[j][count]: # check that it fits 100%
                    # print(j,ref_index,count)
                    # close_index = resi_to_index(closeAA,structure_list[j],structure_list,align)
                    # if close_index != None:
                    print(j, closeAA)
                    if bigger_AA(seq[closeAA],align[0][closeAA]):
                        flag = False     
                if flag:
                    residue = resi
                    resi_list.append(index_to_resi(ref_index,seq,atoms_list[k]))
                else:
                    resi_list.append("-")
            else:
                resi_list.append("-")
        else:
            resi_list.append("-")
    # except:
    #     resi_list.append("-")
    return residue, resi_list


# def discard_surface_residues(hotspot,structure):
#     import findsurfaceatoms
#     exposed_residues = findsurfaceatoms.findSurfaceResidues(selection=structure)
#     exposed_set = set()
#     for resi in exposed_residues:
#         exposed_set.add(resi[1])
#     new_hotspot = {}
#     for k,v in hotspot.items():
#         if k not in exposed_set:
#             new_hotspot[k] = v
#     return new_hotspot, exposed_set



# def get_close_aa(close_AAs,modelatoms,kd,atom,resi):
#     tmp_set = set([int(modelatoms[x].resi) for x in kd.query_ball_point(atom.coord,4)])
#     tmp_set.discard(resi)
#     close_AAs.update(tmp_set)
#     return close_AAs


# def get_close_aa_list(structure_list,align):
#     closeAA_list = []
#     for i, seq in enumerate(align):
#         closeAA_list.append([])
#         for ele in seq:
#             closeAA_list[i].append(ele)

#     for i, structure in enumerate(structure_list):
#         model = cmd.get_model(structure+" and (not name CA and not name N and not name C and not name O)")
#         kd = cKDTree([atom.coord for atom in model.atom])
#         resi = 0
#         close_AAs = set()
#         for atom in model.atom:
#             if resi != int(atom.resi):
#                 if resi != 0:
#                     # if resi_to_index(resi,structure,structure_list,align) != None:
#                     closeAA_list[i][resi_to_index(resi,structure,structure_list,align)] = close_AAs
#                     # else:
#                         # pass
#                     close_AAs = set()
#                 resi = int(atom.resi)
#                 close_AAs = get_close_aa(close_AAs,model.atom,kd,atom,resi,structure,structure_list,align)
#             else:
#                 close_AAs = get_close_aa(close_AAs,model.atom,kd,atom,resi,structure,structure_list,align)
#         # if resi_to_index(resi,structure,structure_list,align) != None:
#         closeAA_list[i][resi_to_index(resi,structure,structure_list,align)] = close_AAs
#         # else:
#             # pass
#     return closeAA_list


def run(structure_list,alignment_file_name,structure_list_chainA,score_list):
    print("Finding hotspots...")
    align = AlignIO.read(alignment_file_name,"clustal")


    atoms_list = []
    for structure in structure_list:
        atoms_list.append(cmd.get_model(structure+" and name CA").atom)

    # print("start")
    # closeAA_list = get_close_aa_list(structure_list,align)
    # print(closeAA_list)
    closeAA_list, hotspot_involved = get_new_close_AA(structure_list,align,score_list,structure_list_chainA,atoms_list)
    # print("outer",len(closeAA_list))
    # for e in closeAA_list:
    #     print("inner",len(e))
    # print("hej")
    # print(closeAA_list)
    # hotspot_list = []
    # exposed_list = []
    # for i, structure in enumerate(structure_list):
    hotspot = {}
    count = 0
    for atom in atoms_list[0]:
        if int(atom.resi) in hotspot_involved:
            k,v = add_hotspot(closeAA_list,structure_list,align,atoms_list,atom.resn,int(atom.resi),atom.coord,count)
            if k != None:
                hotspot[k] = v
            count += 1
    # print(hotspot)
    # if discard_exposed:
    #     hotspot, exposed_set = discard_surface_residues(hotspot,structure)
    #     exposed_list.append(exposed_set)
    # else:
    #     exposed_set = None
    # hotspot_list.append(hotspot)
    return hotspot


def print_hotspot(hotspot,structure_list):
    model_list = []
    for structure in structure_list:
        model_list.append(cmd.get_model(structure+" and name CA").atom)
    # for i, hotspot in enumerate(hotspot_list):
    print("\tPrinting possible single mutations in "+structure_list[0])
    for k,v in hotspot.items():
        for j, resi in enumerate(v):
            if resi != "-":
                print("\t\t"+model_list[0][k-1].resn+model_list[0][k-1].resi+" -> "+model_list[j+1][resi-1].resn+" as structure: "+structure_list[j+1])
    # if print_hospots_from_structure == "ref_structure":
    #     break


from statistics import median


# for i, structure in enumerate(structure_list):
#     model = cmd.get_model(structure+" and (not name CA and not name N and not name C and not name O)")
#     kd = cKDTree([atom.coord for atom in model.atom])
#     resi = 0
#     close_AAs = set()
#     for atom in model.atom:
#         if resi != int(atom.resi):
#             if resi != 0:
#                 # if resi_to_index(resi,structure,structure_list,align) != None:
#                 closeAA_list[i][resi_to_index(resi,structure,structure_list,align)] = close_AAs
#                 # else:
#                     # pass
#                 close_AAs = set()
#             resi = int(atom.resi)
#             close_AAs = get_close_aa(close_AAs,model.atom,kd,atom,resi)
#         else:
#             close_AAs = get_close_aa(close_AAs,model.atom,kd,atom,resi)
#     # if resi_to_index(resi,structure,structure_list,align) != None:
#     closeAA_list[i][resi_to_index(resi,structure,structure_list,align)] = close_AAs
#     # else:
#         # pass

# def get_close_aa(close_AAs,modelatoms,kd,atom,resi):
#     tmp_set = set([int(modelatoms[x].resi) for x in kd.query_ball_point(atom.coord,4)])
#     tmp_set.discard(resi)
#     return close_AAs.union(tmp_set)

# def resi_set_to_index(resi_set,ref_structure,structure_list,align):
#     output = set()
#     for ele in resi_set:
#         output.add(resi_to_index(ele,ref_structure,structure_list,align))
#         # output.add(ele)
#     # print("getClose",close_AAs)
#     return output


def update_closeAA_list(close_AA_list, tmp_close_AAs, count, residue, resi_list, index, align_seq, atomsCA):
    close_AAs = set()
    for ele in tmp_close_AAs:
        close_AAs.add(resi_to_index(ele,align_seq, atomsCA))
    try:
        while len(resi_list[:resi_list.index(residue)]) > count:
            # close_AA_list[j][count] = close_AAs
            close_AA_list[index].append(close_AAs)
            count += 1
            close_AAs = set()
    except:
        pass
    return close_AA_list, close_AAs, count

def get_new_close_AA(structure_list,align,score_list,structure_list_chainA,atoms_list):

    # Get all residues of upper median and the ones within 4Å of these in the reference structure
    # models = []
    # kds = []
    # for i, structure in enumerate(structure_list_chainA):    # Get models
    #     model = cmd.get_model(structure+" and (not name CA and not name N and not name C and not name O)")
    #     models.append(model)
    #     kds.append(cKDTree([atom.coord for atom in model.atom]))
    # atoms = cmd.get_model(structure_list_chainA[0]+" and (not name CA and not name N and not name C and not name O)").atom
    # CAatoms = cmd.get_model(structure_list_chainA[0]+" and name CA").atom
    # kd = cKDTree([atom.coord for atom in atoms])
    CAatoms = atoms_list[0]
    kdCA = cKDTree([atom.coord for atom in CAatoms])
    ref_structure = structure_list[0]
    hotspot_involved_set = set()
    med = median(score_list[0])
    for i, score in enumerate(score_list[0]):    # Get upper median
        if score > med:
            hotspot_involved_set.add(int(CAatoms[i].resi))
    # print(hotspot_involved_set)
    tmp = set()
    # kage = ""
    for atom in CAatoms:     # Add close AA to hospot_involved_set
        residue = int(atom.resi)
        # if True:
        if residue in hotspot_involved_set:
            # print(residue)
            # # print(set([int(atoms[x].resi) for x in kd.query_ball_point(atom.coord,4)]))
            # print(set([int(CAatoms[x].resi) for x in kdCA.query_ball_point(atom.coord,6)]))
            # if len(set([int(CAatoms[x].resi) for x in kdCA.query_ball_point(atom.coord,6)])) > 7:
            #     print(residue)
            #     kage += f" or resi {residue}"
            tmp.update(set([int(CAatoms[x].resi) for x in kdCA.query_ball_point(atom.coord,6)]))
            # print(residue,set([int(CAatoms[x].resi) for x in kdCA.query_ball_point(atom.coord,6)]))
            # print(tmp)
    # print(tmp)
    # print(kage)
    # print(len(tmp))
    # for residue in hotspot_involved_set:     # Add residues within 4Å to upper median residues
    #     # atoms = models[0].atom
    #     ref_atom = atoms[int(residue-1)]
    #     tmp.union(set([int(atoms[x].resi) for x in kd.query_ball_point(ref_atom.coord,4)]))
    # print(tmp)
    hotspot_involved_set.update(tmp)
    # print(len(hotspot_involved_set))
    # print(hotspot_involved_set)
    hotspot_involved_index = [] 
    # hotspot_involved_resi = []      # Change hotspot_involved to list with alignment index
    for residue in hotspot_involved_set:
        hotspot_involved_index.append(resi_to_index(residue,align[0],CAatoms))
        # hotspot_involved_resi.append(residue)
    # hotspot_involved_set = set(hotspot_involved)
    # print(hotspot_involved)
    # # ^^^^^ FINE!
    # print(len(hotspot_involved))

    # Make close_AA_list -> [[set(), set()]]. set is close AA to given hotspot_involved. Inner list is len of hotspot_involved. Outer_list is len of structure_list - ref_structure
    close_AA_list = []
    for j, structure in enumerate(structure_list[1:]):
        # print(structure)
        resi_list = [0]
        for index in hotspot_involved_index:
            resi_list.append(index_to_resi(index,align[j+1],atoms_list[j+1]))
        resi_set = set(resi_list)
        # print(structure,resi_set)
        close_AA_list.append([])
        atoms = cmd.get_model(structure_list_chainA[j+1]+" and (not name CA and not name N and not name C and not name O)").atom
        # print("hejhehj",structure_list_chainA[j+1])
        kd = cKDTree([atom.coord for atom in atoms])
        residue = 0
        # ref_index = 0
        # kage = []
        close_AAs = set()
        count = 0
        for atom in atoms:
            if residue != int(atom.resi):
                close_AA_list, close_AAs, count = update_closeAA_list(close_AA_list, close_AAs, count, residue, resi_list, j, align[j+1], atoms_list[j+1])
                residue = int(atom.resi)
            if residue in resi_set:
                tmp_set = set([int(atoms[x].resi) for x in kd.query_ball_point(atom.coord,4)])
                tmp_set.discard(residue)
                close_AAs.update(tmp_set)
        close_AA_list, close_AAs, count = update_closeAA_list(close_AA_list, close_AAs, count, residue, resi_list, j, align[j+1], atoms_list[j+1])
        # print(j+1,close_AA_list[j])
        #     if residue != int(atom.resi):
        #         if tmp != "":
        #             if residue + 1 != int(atom.resi) and :
        #                 # if tmp != "":
        #                     close_AA_list[j][count] = "G"
        #                     # count += 1
        #                     kage.append(atom.resi)
        #             else:
        #                 # if tmp != "":
        #                     close_AA_list[j][count] = resi_set_to_index(tmp,structure,structure_list,align)
        #                     # count += 1
        #                     kage.append(atom.resi)
        #             count += 1
        #         tmp = set()
        #         residue = int(atom.resi)
        #         if residue in resi_set:
        #             # if tmp != "":
        #             #     close_AA_list[j][count] = resi_set_to_index(tmp,structure,structure_list,align)
        #             #     count += 1
        #             #     kage.append(atom.resi)
        #             # tmp = set()
        #             tmp = get_close_aa(tmp,atoms,kd,atom,residue)
        #         # else:
        #         #     print(atom.resi)
        #         # if residue != 0:
        #         #     residue = int(atom.resi)
        #         #     print(j,count,residue)
        #         #     close_AA_list[j][count] = resi_set_to_index(tmp,structure,structure_list,align)
        #         #     count += 1
        #         #     tmp = set()
                
        #         # # print("NEW",residue,structure)
        #         # if residue in resi_set:
        #         #     # print("NEW",tmp)
        #         #     tmp = get_close_aa(tmp,atoms,kd,atom,residue)
        #         # print("je",atom.resi)
        #         # if residue in resi_set:  #NU
        #         #     print("jjje",atom.resi)
        #         #     if residue != 0:      #OLD
                        
        #         #         close_AA_list[j][count] = tmp
        #         #         count += 1
        #         #         tmp = set()
        #         #     residue = int(atom.resi)    #UPDATEresi
        #         #     tmp = get_close_aa(tmp,atoms,kd,atom,residue,structure,structure_list,align) #

                    
        #         # ref_index = resi_to_index(residue,structure_list[j],structure_list,align)
        #     elif residue in resi_set:  #NU
        #         # print("SAME",residue)
        #         # print("SAME",tmp)
        #         tmp = get_close_aa(tmp,atoms,kd,atom,residue)
        #     # if ref_index in hotspot_involved_set:
        # close_AA_list[j][count] = resi_set_to_index(tmp,structure,structure_list,align)
        # kage.append(atom.resi)
        # for i, ele in enumerate(resi_list):
        #     if ele != kage[i]:
        #         print("HALLLOOOO",ele,kage[i])
        # for i, index in enumerate(hotspot_involved):
        #     tmp = set([int(atoms[x].resi) for x in kd.query_ball_point(ref_atom.coord,4)])
        #     tmp.discard(index_to_resi(index,ref_structure,structure_list,align))
        #     close_AA_list[j][i] = tmp
    return close_AA_list, hotspot_involved_set


# def get_close_aa(close_AAs,modelatoms,kd,atom,resi):
#     tmp_set = set([int(modelatoms[x].resi) for x in kd.query_ball_point(atom.coord,4)])
#     tmp_set.discard(resi)
#     return close_AAs.union(tmp_set)


# for i, structure in enumerate(structure_list):
#     model = cmd.get_model(structure+" and (not name CA and not name N and not name C and not name O)")
#     kd = cKDTree([atom.coord for atom in model.atom])
# #     resi = 0
# #     close_AAs = set()
#     for atom in model.atom:
#         if resi != int(atom.resi):
#             if resi != 0:
#                 # if resi_to_index(resi,structure,structure_list,align) != None:
#                 closeAA_list[i][resi_to_index(resi,structure,structure_list,align)] = close_AAs
#                 # else:
#                     # pass
#                 close_AAs = set()
#             resi = int(atom.resi)
#             close_AAs = get_close_aa(close_AAs,model.atom,kd,atom,resi)
#         else:
#             close_AAs = get_close_aa(close_AAs,model.atom,kd,atom,resi)
#     # if resi_to_index(resi,structure,structure_list,align) != None:
#     closeAA_list[i][resi_to_index(resi,structure,structure_list,align)] = close_AAs
#     # else:
#         # pass
    
    
    # resi = 0
    # for atom in model.atom:
    #     if resi != int(atom.resi):
    #         if resi != 0:
    #             hotspot_involved

    
#    model_kd = dict()  
#         for structure in structure_list_entire_chainA:
#             model = cmd.get_model(structure+" and name CA")
#             model_kd[model] = cKDTree([atom.coord for atom in model.atom])