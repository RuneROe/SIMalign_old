from pymol import cmd
import numpy as np
from scipy.spatial import cKDTree
import re
import sys


def aa_to_blosom(aa1,aa2):
    """
    DESCRIPTION

    Input: 2 amino acid residues (UPPER case tree-letter-code)
    output: blosom62 score

    table from https://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
    """
    blosom62 = [
    [4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0,-2,-1,0,-4],
    [-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3,-1,0,-1,-4],
    [-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3,3,0,-1,-4],
    [-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3,4,1,-1,-4],
    [0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4],
    [-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2,0,3,-1,-4],
    [-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2,1,4,-1,-4],
    [0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3,-1,-2,-1,-4],
    [-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3,0,0,-1,-4],
    [-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3,-3,-3,-1,-4],
    [-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1,-4,-3,-1,-4],
    [-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2,0,1,-1,-4],
    [-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1,-3,-1,-1,-4],
    [-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1,-3,-3,-1,-4],
    [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2,-2,-1,-2,-4],
    [1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2,0,0,0,-4],
    [0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0,-1,-1,0,-4],
    [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3,-4,-3,-2,-4],
    [-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1,-3,-2,-1,-4],
    [0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4,-3,-2,-1,-4],
    [-2,-1,3,4,-3,0,1,-1,0,-3,-4,0,-3,-3,-2,0,-1,-4,-3,-3,4,1,-1,-4],
    [-1,0,0,1,-3,3,4,-2,0,-3,-3,1,-1,-3,-1,0,-1,-3,-2,-2,1,4,-1,-4],
    [0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2,0,0,-2,-1,-1,-1,-1,-1,-4],
    [-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,1]]
    aa_id_list = ["ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","B","Z","X","-"]
    aa_to_id = dict()
    for i, aa in enumerate(aa_id_list):
        aa_to_id[aa] = i
    id1 = aa_to_id[aa1]
    id2 = aa_to_id[aa2]
    return blosom62[id1][id2]



def select_by_score(score_list, modelatoms):
    """
    DESCRIPTION

    Takes list of scores and make a list of strings used for selection for the residues with score above procentive tresshold

    DEPENDENCIES

    import numpy as np
    """

    tresshold = np.median(np.array(list(filter(lambda num: num != 0, score_list))))
    out_string = ""
    chain = ""
    for j, score in enumerate(score_list):
        if score > tresshold:
            if chain == "":
                chain = modelatoms[j].chain
                out_string = f" and ((chain {chain} and (resi {j+1}"
            elif chain != modelatoms[j].chain:
                chain = modelatoms[j].chain
                out_string += f")) or (chain {chain} and (resi {j+1}"
            else:
                out_string += f" or resi {j+1}"
    return out_string+")))"
    

def downloading_files(ref_structure,files):
    try:
        cmd.load(ref_structure)
        for file in files:
            if file != ref_structure:
                if len(file.split(".")) == 1:
                    print("\t\tFetch:",file)
                    cmd.fetch(file)
                else:
                    print("\t\tLoad:",file)
                    cmd.load(file)
        structure_list_entire = cmd.get_object_list()
        
        return structure_list_entire[0], structure_list_entire
    except:
        print("\n- - - - - - - - - - - - - - - - - -\nImport ERROR")
        print("\nCouldn't import one or more of the files in pymol\n- - - - - - - - - - - - - - - - - -\n")
        sys.exit(1)
    
from Bio import AlignIO



def get_align(alignment_file_name,structure_list):
    alignIO = AlignIO.read(alignment_file_name,"clustal")
    align = []
    resi = []
    for seq in alignIO:
        resi.append(0)
    modelsatoms = [cmd.get_model(structure+" and chain A and not HETATM and name CA").atom for structure in structure_list]
    for i in range(len(alignIO[0].seq)):
        tmp = {}
        for j, seq in enumerate(alignIO):
            if seq[i] != "-":
                try:
                    tmp[structure_list[j]] = modelsatoms[j][resi[j]].index
                except IndexError:
                    pass
                resi[j] += 1
        align.append(tmp)
    return align


def update_alignment(align):
    new_align = []
    for pos in align:
        tmp = []
        for k,v in pos.items():
            tmp.append((k,v))
        new_align.append(tmp)
    cmd.set_raw_alignment("aln",new_align)

def process_dicts(dict_list, namelist):
    dictlen = len(dict_list)
    for i in range(len(namelist)):
        for j in range(i + 1, len(namelist)):
            for dict1 in dict_list:
                for dict2 in dict_list:
                    if dict1 != dict2:
                        try:
                            if dict1[namelist[i]] < dict2[namelist[i]] and dict1[namelist[j]] > dict2[namelist[j]]:
                                if abs(dict2[namelist[i]]-dict2[namelist[j]]) > dictlen:
                                    del dict2[namelist[j]]
                                else:
                                    del dict1[namelist[j]]
                        except:
                            pass
    return dict_list

def average_coordinate(list_of_coordinates):
    n = len(list_of_coordinates)
    for i, x in enumerate(list_of_coordinates):
        if i == 0:
            average = np.array(x)
        else:
            average += np.array(x)
    average = average/n
    return average


def SIMalign(ref_structure, structure_list_entire, iterations, tresshold_aa, max_dist, max_initial_rmsd):
    n_homologous_list = len(structure_list_entire) - 1
    break_flag = False
    to_outfile = [""]
    selection = []

    for j in range(iterations):

        to_outfile.append(f"\tIteration {j+1}\n\tstructure\tRMSD\tresidues aligned\n")
        tmp_out = ""
        
        # Super of all structures to ref_structure
        for i, structure in enumerate(structure_list_entire):
            if i != 0:
                print(f"\tsuperimposing",structure,"towards",ref_structure)
                if selection == []:
                    super = cmd.super(target=f"{ref_structure} and name CA and not HETATM and chain A", mobile=f"{structure} and name CA and not HETATM and chain A")
                    if super[0] > max_initial_rmsd:
                        print(f"\tStructure {structure} was deleted because the RMSD to {ref_structure} was above 5Å: {super[0]}Å")
                        cmd.delete(structure)
                else:
                    super = cmd.super(target=f"{ref_structure} and name CA and not HETATM{selection[0]}", mobile=f"{structure} and name CA and not HETATM{selection[i]}")
                tmp_out += f"\t{structure_list_entire[i]}\t{round(super[0],3)}\t{super[1]}\n"
        to_outfile.append(tmp_out)
        selection = []
        structure_list_entire = cmd.get_object_list()

        # Get models and cKDtree
        model_kd = dict()  
        for structure in structure_list_entire:
            model = cmd.get_model(structure+" and name CA and not HETATM and chain A")
            model_kd[model] = cKDTree([atom.coord for atom in model.atom])
        
        score_list = []
        
        align = []

        # LOOP through models
        for i, model in enumerate(model_kd.keys()):

            score = []
            ref_kd = model_kd[model]

            # LOOP through ca atoms to find scores of one model
            for ref_atom in model.atom:
                
                s = 0
                ref_resn = ref_atom.resn
                ref_coord = ref_atom.coord
                max_score = aa_to_blosom(ref_resn,ref_resn) + 4

                tmp_coordinates = [ref_coord]
                # Finding closest AA from other models to ref AA
                for x, (m, kd) in enumerate(model_kd.items()):
                    
                    if m != model:
                        closest_pair = kd.query(ref_coord)
                        atom = m.atom[closest_pair[1]]
                        if closest_pair[0] == ref_kd.query(atom.coord)[0] and closest_pair[0] <= max_dist: 
                            s += ((aa_to_blosom(ref_resn,atom.resn) + 4)/(n_homologous_list*max_score))
                            tmp_coordinates.append(atom.coord)
                
                # Updating seq alignment so it is based on structureal information
                tmp_center = average_coordinate(tmp_coordinates)
                tmp = {}
                for x, (m, kd) in enumerate(model_kd.items()):
                    closest_pair = kd.query(tmp_center)
                    atom = m.atom[closest_pair[1]]
                    if closest_pair[0] <= max_dist: 
                        tmp[structure_list_entire[x]] = atom.index
                score.append(s)
            
                if tmp not in align:
                    flag = True
                    for ele in align:
                        for structure in structure_list_entire:
                            if structure in ele and structure in tmp:
                                if tmp[structure] == ele[structure]:
                                    flag = False
                    if flag:
                        if align == []:
                            align.append(tmp)
                        else:
                            for structure in structure_list_entire:
                                if structure in tmp.keys():
                                    key = structure
                                    break
                            for x in range(len(align)+1):

                                try:
                                    if x == len(align):
                                        align.append(tmp)
                                        break
                                    elif tmp[key] < align[x][key]:
                                        align = align[:x]+[tmp]+align[x:]
                                        break
                                except:
                                    pass

           
            score_list.append(score)
            tmp_selection = select_by_score(score, model.atom)
            # Checking that it is not below tresshold_aa
            if len(re.findall(r"\s\d", tmp_selection)) < tresshold_aa:
                print("\t"+structure_list_entire[i], len(re.findall(r"\s\d", tmp_selection)))
                break_flag = True
            selection.append(tmp_selection)
        
        if break_flag:
            print(f"\tBreaked after {j} iteration(s) of superexposion. \nTry to change the parameter tresshold_aa if a higher number of iterations are wanted.")
            break

        align = process_dicts(align,structure_list_entire)

        print(to_outfile[-2]+to_outfile[-1][:-1])
        if to_outfile[-1] == to_outfile[-3]:
            break_flag = True
            print(f"\tBreaked after {j+1} iteration(s) of superimposion.")
            break
        
    with open("log.txt","w") as outfile:
        outfile.write("".join(to_outfile))
    update_alignment(align)
    if break_flag == False:
        print(f"\tCompleted {iterations} iteration(s) of superimposion.")
    return score_list, selection, structure_list_entire

def run(ref_structure, files, iterations, tresshold_aa, max_dist, alignment_file_name, max_initial_rmsd):
    """
    DESCRIPTION

    Main function that takes files containing structures and return a .pse file where the structures are aligned and colored after similarity.

    DEPENDENCIES

    from pymol import cmd
    import numpy as np
    from scipy.spatial import cKDTree
    import re
    import sys
    select_by_score
    color_pymol
    color_by_number
    """
    cmd.reinitialize()
    print("Loading structures to pymol...")
# Importing files - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - -
    ref_structure, structure_list_entire = downloading_files(ref_structure,files)
    # cmd.alignto(ref_structure)

# Basic Pymol stuff - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - -

    cmd.remove("hydrogens")


# LOOP start
    print("Running SIMalign...")
    score_list, core_selection, structure_list_entire = SIMalign(ref_structure, structure_list_entire, iterations, tresshold_aa, max_dist, max_initial_rmsd)
    cmd.save(alignment_file_name, selection="aln")

    

    return len(cmd.get_model(ref_structure).atom), score_list, structure_list_entire, core_selection