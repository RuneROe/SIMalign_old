
#TO DO
#- Lav bedre output fra run function
#- Lave visualisering i ipynb med nyt py doc
#- Lav et print der opdaterer sig selv, så man kan have et fedt download.



from pymol import cmd
import numpy as np
from scipy.spatial import cKDTree
import re
import sys


#From https://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
def aa_to_blosom(aa1,aa2):
    """
    DESCRIPTION

    Input: 2 amino acid residues (UPPER case tree-letter-code)
    output: blosom62 score
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
        # if remove_chain_duplicate:
        #     ref_structure = ref_structure.split(".")[0]+"__CA"
        # else:
        #     ref_structure = ref_structure.split(".")[0]+"_CA"
        for file in files:
            if file != ref_structure:
                cmd.load(file)
        structure_list_entire = cmd.get_object_list()
        # if remove_chain_duplicate:
        #     cmd.hide("everything", "all and not HETATM")
        #     for i, structure in enumerate(structure_list_entire):
        #         tmp_chains = cmd.get_fastastr(structure).split(">")[1:]
        #         chains = dict()
        #         for c in tmp_chains:
        #             c = c.split("\n")
        #             value = "".join(c[1:])
        #             if value not in chains.values():
        #                 chains[c[0]] = value
        #         cmd.select(f"{structure}_", "none")
        #         for k in chains.keys():
        #             chain_id = k.split("_")[-1]
        #             cmd.show("cartoon", f"{structure} and chain {chain_id}")
        #             cmd.select(f"{structure}_", f"{structure}_ or (chain {chain_id} and {structure})")
        #         structure_list_entire[i] = f"{structure}_"
        
        return structure_list_entire
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
        print(len("".join([str(x) for x in seq.seq.split("-")])))
    modelsatoms = [cmd.get_model(structure+" and chain A and not HETATM and name CA").atom for structure in structure_list]
    print(cmd.get_object_list())
    for ele in modelsatoms:
        print(len(ele))
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


# def update_seq_alignment(align,tmp,i,n_homologous_list,structure_list,ref_structure,ref_index):
#     align_set = set([set(x) for x in align])
#     if set(tmp) not in align_set and i != n_homologous_list:
#         flag = True
#         for i,tmp_set in enumerate(align):
#             if (ref_structure, ref_index) in tmp_set:
#                 break
#         for ele in structure_list[:i]:
#             #if not gab
#                 flag = False
#                 break
#         if flag:
#             count = 0
#             for ele in structure_list[i+1:]:
#                 count += 1
#                 for wt in structure_list[i+count:]:
#                     pass #check if  they are in align
#                         flag = False
#                         break
#             # update align
#     return align

def update_alignment(align):
    new_align = []
    for pos in align:
        tmp = []
        for k,v in pos.items():
            tmp.append((k,v))
        new_align.append(tmp)
    cmd.set_raw_alignment("aln",new_align)
        
def SIMalign(ref_structure, structure_list_entire, iterations, tresshold_aa, max_dist, alignment_file_name):
    n_homologous_list = len(structure_list_entire) - 1
    # align_structure_list = []
    # for i, structure in enumerate(structure_list_entire):
    #     align_structure_list.append(structure+"_align")
    break_flag = False
    to_outfile = [""]
    for j in range(iterations):
        align = get_align(alignment_file_name,structure_list_entire)
        # Get models and cKDtree
        model_kd = dict()  
        for structure in structure_list_entire:
            model = cmd.get_model(structure+" and name CA and not HETATM and chain A")
            model_kd[model] = cKDTree([atom.coord for atom in model.atom])
        
        score_list = []
        selection = []

        # LOOP through models
        for i, model in enumerate(model_kd.keys()):

            score = []
            ref_kd = model_kd[model]

            # LOOP through ca atoms to find scores of one model
            for ref_atom in model.atom:
                tmp = {structure_list_entire[i]: ref_atom.index}
                s = 0
                ref_resn = ref_atom.resn
                ref_coord = ref_atom.coord
                max_score = aa_to_blosom(ref_resn,ref_resn) + 4

                # Finding closest AA from other models to ref AA
                for x, (m, kd) in enumerate(model_kd.items()):
                    
                    if m != model:
                        closest_pair = kd.query(ref_coord)
                        atom = m.atom[closest_pair[1]]
                        if closest_pair[0] == ref_kd.query(atom.coord)[0] and closest_pair[0] <= max_dist: 
                            s += ((aa_to_blosom(ref_resn,atom.resn) + 4)/(n_homologous_list*max_score))
                            tmp[structure_list_entire[x]] = atom.index

                #changing alignment
                if tmp not in align:
                    for ele in align:
                        if structure_list_entire[i] in ele:
                            align.remove(ele)
                            align.append(tmp)
                    for t in tmp.keys():
                        if t != structure_list_entire[i]:
                            for ele2 in align:
                                if ele2 != tmp and t in ele2:
                                    align.remove(ele2)


                score.append(s)
            score_list.append(score)

            tmp_selection = select_by_score(score, model.atom)
            # Checking that it is not below tresshold_aa
            if len(re.findall(r"\s\d", tmp_selection)) < tresshold_aa:
                print(structure_list_entire[i], len(re.findall(r"\s\d", tmp_selection)))
                break_flag = True
            selection.append(tmp_selection)
        if break_flag:
            print(f"Breaked after {j} iteration(s) of superexposion. \nTry to change the parameter tresshold_aa if a higher number of iterations are wanted.")
            break
        to_outfile.append(f"Iteration {j+1}\nstructure\tRMSD\tatoms aligned\n")
        tmp_out = ""
        for i, structure in enumerate(structure_list_entire):
            # cmd.select(align_structure_list[i], structure+selection[i])
            if i != 0:
                # super = cmd.super(target=align_structure_list[0], mobile=align_structure_list[i])
                super = cmd.super(target=f"{ref_structure} and name CA and not HETATM and chain A{selection[0]}", mobile=f"{structure} and name CA and not HETATM and chain A{selection[i]}")
                tmp_out += f"{structure_list_entire[i]}\t{round(super[0],3)}\t{super[1]}\n"
        to_outfile.append(tmp_out)
        print(to_outfile[-2]+to_outfile[-1])
        if to_outfile[-1] == to_outfile[-3]:
            break_flag = True
            print(f"Breaked after {j+1} iteration(s) of superexposion.")
            break
        
    with open("log.txt","w") as outfile:
        outfile.write("".join(to_outfile))
    update_alignment(align)
    if break_flag == False:
        print(f"Completed {iterations} iteration(s) of superexposion.")
    for i, structure in enumerate(structure_list_entire):
        cmd.select(f"{structure}_core",f"{structure} and not HETATM and chain A{selection[i]}")
    return score_list

def run(ref_structure, files, iterations, tresshold_aa, max_dist, alignment_file_name):
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
    structure_list_entire = downloading_files(ref_structure,files)
    ref_structure = ref_structure.split(".")[0]
    

# Basic Pymol stuff - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - -

    cmd.remove("hydrogens")
    cmd.alignto(ref_structure+" and chain A", object="aln")
    cmd.save(alignment_file_name, selection="aln")
    print("saved")

# LOOP start
    print("Running SIMalign...")
    score_list = SIMalign(ref_structure, structure_list_entire, iterations, tresshold_aa, max_dist, alignment_file_name)
    cmd.save(alignment_file_name, selection="aln")


# Some pymol stuff - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - -
    cmd.hide("cgo", "aln")
    cmd.set("seq_view_label_mode", "1")
    cmd.set("antialias", "4")
    cmd.set("ray_trace_mode", "1")

    # cmd.save(outfilename)
    

    return len(cmd.get_model(ref_structure).atom), score_list, structure_list_entire





