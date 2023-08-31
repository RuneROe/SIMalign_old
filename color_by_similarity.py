
#TO DO
#- Lav bedre output fra run function
#- rydde op i py doc
#- Lave visualisering i ipynb med nyt py doc
#- Lav et print der opdaterer sig selv, så man kan have et fedt download.



from pymol import cmd
import numpy as np
import os
import sys
import time
from scipy.spatial import cKDTree
import re

start_time = time.time()


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



def color_by_number(number):
    """
    DESCRIPTION

    Rainbow from red to white
    number 0: red
    number 1: white

    DEPENDENCIES

    import numpy as np
    """
    return [0.8+(number/5),number,number]


def color_pymol(model_CA, structure_entire, score_list):
    """
    DESCRIPTION



    DEPENDENCIES

    color_by_number
    from pymol import cmd
    """
    # model = cmd.get_model(structure_entire)
    for i, atom in enumerate(model_CA.atom):
        print(structure_entire, atom.resi, atom.resn, atom.chain, score_list[i], color_by_number(score_list[i]))
        cmd.set_color(f"{atom}color", color_by_number(score_list[i]))
        cmd.color(f"{atom}color", f"resi {atom.resi} and chain {atom.chain} and {structure_entire}")

    # cmd.set_color(f"color{ref_atom}", color_by_number(number))
    # cmd.color(f"color{ref_atom}", f"resi {ref_atom.resi} and {ref_structure}")

def select_by_score(score_list, modelatoms):
    """
    DESCRIPTION

    Takes list of scores and make a list of strings used for selection for the residues with score above procentive tresshold

    DEPENDENCIES

    import numpy as np
    """
    # score_list_no_zeroes = np.array(list(filter(lambda num: num != 0, score_list)))
    # max = np.max(score_list)
    tresshold = np.median(np.array(list(filter(lambda num: num != 0, score_list))))
    out_string = ""
    chain = ""
    # tresshold = (max+median)/2
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
    




def run(ref_structure, files, iterations, tresshold_aa, max_dist, remove_chain_duplicate, outfile):
# redundency_tresshold = 0.8 # Removes structure that are more similar than the tresshold to others. 1 means that nothing is removed
# User arguments - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# arguments = sys.argv
# if len(arguments) == 1:
#     path_ref_structure = input("Please type the path for the reference structure: ")
#     path_folder = input("Please type the path for the folder with the homologues structures: ")
# elif len(arguments) == 2:
#     path_ref_structure = arguments[1]
#     path_folder = input("Please type the path for the folder with the homologues structures: ")
# else:
#     path_ref_structure = arguments[1]
#     path_folder = arguments[2]
#     if len(arguments) > 3:
#         try:
#             iterations = int(arguments[3])
#         except:
#             print("\n- - - - - - - - - - - - - - - - - -\nKEY ERROR")
#             print("iterations needs to be typed as an int!")
#             print("\nUSAGE\n\ncolor_by_similarity.py ref_structure.pdb homologues_structure_folder [iterations=3] [tresshold_aa=100]\n- - - - - - - - - - - - - - - - - -\n")
#             sys.exit(1)
#         if len(arguments) > 4:
#             try:
#                 tresshold_aa = float(arguments[4])
#             except:
#                 print("\n- - - - - - - - - - - - - - - - - -\nKEY ERROR")
#                 print("Tresshold aa needs to be typed as a float!")
#                 print("\nUSAGE\n\ncolor_by_similarity.py ref_structure.pdb homologues_structure_folder [iterations=3] [tresshold_aa=100]\n- - - - - - - - - - - - - - - - - -\n")
#                 sys.exit(1)

# path_ref_structure = "ublipase.pdb"
# path_ref_structure = "ref_structure.pdb"
# path_folder = "homologues_plastic"



# Importing files - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #try:
        # cmd.fetch("2gub")
        # cmd.fetch("3a6z")
        # cmd.fetch("2zvd")
        # cmd.fetch("3a70")
        # ref_structure = "ubli"
    cmd.reinitialize
    cmd.load(ref_structure)
    files.remove(ref_structure)
    if remove_chain_duplicate:
        ref_structure = ref_structure.split(".")[0]+"__CA"
    else:
        ref_structure = ref_structure.split(".")[0]+"_CA"
    for file in files:
        cmd.load(file)
    structure_list_entire = cmd.get_object_list()
    if remove_chain_duplicate:
        cmd.hide("everything", "all and not HETATM")
        for i, structure in enumerate(structure_list_entire):
            tmp_chains = cmd.get_fastastr(structure).split(">")[1:]
            chains = dict()
            for c in tmp_chains:
                c = c.split("\n")
                value = "".join(c[1:])
                if value not in chains.values():
                    chains[c[0]] = value
            cmd.select(f"{structure}_", "none")
            for k in chains.keys():
                chain_id = k.split("_")[-1]
                cmd.show("cartoon", f"{structure} and chain {chain_id}")
                cmd.select(f"{structure}_", f"{structure}_ or (chain {chain_id} and {structure})")
            structure_list_entire[i] = f"{structure}_"
    n_homologous_list = len(structure_list_entire) - 1
    # ref_model = cmd.get_model(f"{ref_structure} and name CA and resi 1:439 and chain A")
    # homologous_models = []
    structure_list = []
    align_structure_list = []
    for i, structure in enumerate(structure_list_entire):
        structure_list.append(structure+"_CA")
        align_structure_list.append(structure+"_align")
        cmd.select(structure_list[i], structure+" and name CA and not HETATM")
        # homologous_models.append(cmd.get_model(f"{structure} and name CA"))
    # models = [ref_model] + homologous_models
    # print(models)
    # except:
    #     print("\n- - - - - - - - - - - - - - - - - -\nImport ERROR")
    #     print("\nCouldn't import one of the files in pymol\n- - - - - - - - - - - - - - - - - -\n")
    #     sys.exit(1)









# Basic Pymol stuff - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - -

    cmd.remove("hydrogens")



# Functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -





# def scores_by_conservation(ref_model, homologous_list, max_dist):
#     """
#     DESCRIPTION

#     Scoring ref_model acording to sequence similarity to homologous_list with a maximum of max_dist Å between CA.
#     Score 0: unconserved
#     Score 1: conserved

#     DEPENDENCIES

#     aa_to_blosom
    
#     """

#     n_homologous_list = len(homologous_list)
#     score_list = []
#     homologous_residues = []
#     for k, ref_atom in enumerate(ref_model.atom):
#         score = 0
#         ref_resn = ref_atom.resn
#         ref_index = np.array(ref_atom.coord)
#         max_score = aa_to_blosom(ref_resn,ref_resn) + 4
#         dist_list = []
#         tmp = [ref_atom.resi, ref_atom.resn, ref_atom.chain]
#         for model in homologous_list:
#             dist_resn = [model, float("inf"),"-",""]
#             for atom in model.atom:
#                 atom_index = np.array(atom.coord)
#                 dist = np.sqrt(np.sum((atom_index-ref_index)**2))
                
#                 if dist < dist_resn[1]:
#                     dist_resn = [model, dist, atom.resn, atom.resi]
#             if dist_resn[1] <= max_dist:
#                 score += ((aa_to_blosom(ref_resn,dist_resn[2]) + 4)/(n_homologous_list*max_score))
#                 tmp.append(dist_resn[1:])
#             dist_list.append(dist_resn)
#         print(tmp,score)
#         score_list.append(score)
#         homologous_residues.append(dist_list)
#     return score_list#, homologous_residues


# def color_by_number(number):
#     if number < 0.5:
#         return [0.8+(number*2)/5,number*2,number*2]
#     else:
#         return [1-(number-0.5)*2, 1-(number-0.5)*2, 1-((number-0.5)*2)/5]


# def color_by_number(number, minmax=[0.2,0.8], lum_list=[0.2,0.95]):
#     """
#     DESCRIPTION

#     Rainbow from red to blue
#     number 0: red
#     number 1: blue

#     DEPENDENCIES

#     import numpy as np
#     """
#     min = minmax[0]
#     max = minmax[1]
#     if number < 0.25:
#         RGB = np.array([max, number*4*(max-min)+min, min])
#     elif number < 0.5:
#         RGB = np.array([(1-4*(number-0.25))*(max-min)+min, max, min])
#     elif number < 0.75:
#         RGB = np.array([min, max, (number-0.5)*4*(max-min)+min])
#     else:
#         RGB = np.array([min, (1-4*(number-0.75))*(max-min)+min, max])
#     lum = lum_list[0] + (lum_list[1]-lum_list[0])*number
#     tmp_lum = (np.max(RGB)+np.min(RGB))/2
#     out = []
#     if lum > 0.5:
#         for num in RGB:
#             out.append(1-((1-num)*((1-lum)/(1-tmp_lum))))
#     else:
#         for num in RGB:
#             out.append(num*(lum/tmp_lum))
#     return out

# def color_by_number(number, minmax=[0.2,0.8], lum_list=[0.5,1]):
#     """
#     DESCRIPTION

#     Rainbow from red to green
#     number 0: red
#     number 1: green (white)

#     DEPENDENCIES

#     import numpy as np
#     """
#     min = minmax[0]
#     max = minmax[1]
#     if number < 0.5:
#         RGB = np.array([max, number*2*(max-min)+min, min])
#     else:
#         RGB = np.array([(1-2*(number-0.5))*(max-min)+min, max, min])
#     lum = lum_list[0] + (lum_list[1]-lum_list[0])*number
#     tmp_lum = (np.max(RGB)+np.min(RGB))/2
#     out = []
#     if lum > 0.5:
#         for num in RGB:
#             out.append(1-((1-num)*((1-lum)/(1-tmp_lum))))
#     else:
#         for num in RGB:
#             out.append(num*(lum/tmp_lum))
#     return out





# print("\n\n\n")
# print(select_by_score([[0.99,0.5,0.8],[0.81,0.82,0.7]],0.8))
# print("\n\n\n")

# #super
# for i, structure in enumerate(structure_list[1:]):
#     cmd.super(structure,ref_structure)






#udtræk models + shortest model
# models = []
# shortest_model = [float("inf"),""]
# for structure in structure_list:
#     model = cmd.get_model(structure)
#     models.append(model)
#     model_length = len(model.atom)
#     print("model_length =", model_length)
#     if model_length < shortest_model[0]:
#         shortest_model = [model_length, structure]

#draft alignment
# for struc in structure_list_entire[1:]:
#     cmd.super(struc, ref_structure)
    cmd.alignto(ref_structure, object="aln")
    # cmd.alignto(shortest_model[1], object="aln")

    print(structure_list)
    #LOOP start
    break_flag = False
    # outfile = open("log.txt","w")
    to_outfile = [f"structure\tRMSD\tatoms aligned\n"]
    # outfile.write(f"structure\tRMSD\tatoms aligned\n")
    for j in range(iterations):
        # models = [] ### Gammel
        model_kd = dict()  ### NY
        for structure in structure_list:
            model = cmd.get_model(structure)
            # models.append(model)  ### Gammel


# Your list of coordinates
# coordinates = [(x1, y1), (x2, y2), ...]  # Replace with your coordinates

# Create a KDTree from your coordinates
# kdtree = cKDTree(coordinates)

# The point to which you want to find the closest point
# target_point = (target_x, target_y)  # Replace with your target point

# Query the KDTree to find the closest point
# closest_point_idx = kdtree.query(target_point)[1]
# closest_point = coordinates[closest_point_idx]


### NYT
            model_kd[model] = cKDTree([atom.coord for atom in model.atom])
        score_list = []
        selection = []
        for i, model in enumerate(model_kd.keys()):
            # homologous_list = models[:i]+models[i+1:]
            score = []
            homologous_residues = []
            for k, ref_atom in enumerate(model.atom):
                s = 0
                ref_resn = ref_atom.resn
                ref_coord = ref_atom.coord
                max_score = aa_to_blosom(ref_resn,ref_resn) + 4
                dist_list = []
                tmp = [ref_atom.resi, ref_atom.resn, ref_atom.chain]
                for m, kd in model_kd.items():
                    if m != model:
                        closest_point_idx = kd.query(ref_coord)[1]
                        atom = m.atom[closest_point_idx]
                        dist_resn = [m, np.sqrt(np.sum((np.array(ref_coord)-np.array(atom.coord))**2)), atom.resn, atom.resi]
                        # print(np.sqrt(np.sum((np.array(ref_coord)-np.array(atom.coord))**2)))
                        if dist_resn[1] <= max_dist:
                            s += ((aa_to_blosom(ref_resn,dist_resn[2]) + 4)/(n_homologous_list*max_score))
                    #     dist_list.append(dist_resn)
                    #     dist_list.append([])
                    # for model in homologous_list:
                    #     dist_resn = [model, float("inf"),"-",""]
                    #     for atom in model.atom:
                    #         atom_index = np.array(atom.coord)
                    #         dist = np.sqrt(np.sum((atom_index-ref_index)**2))
                            
                    #         if dist < dist_resn[1]:
                    #             dist_resn = [model, dist, atom.resn, atom.resi]
                    #     if dist_resn[1] <= max_dist:
                    #         score += ((aa_to_blosom(ref_resn,dist_resn[2]) + 4)/(n_homologous_list*max_score))
                            tmp.append(dist_resn[1:])
                        dist_list.append(dist_resn)
                print(tmp,s)
                score.append(s)
                homologous_residues.append(dist_list)
            score_list.append(score)

            tmp_selection = select_by_score(score, model.atom)
            if len(re.findall(r"\s\d", tmp_selection)) < tresshold_aa:
                break_flag = True
            selection.append(tmp_selection)
        if break_flag:
            print(f"Colored after {j} iteration(s) of superexposion. \nTry to change the parameter tresshold_aa if a higher number of iterations are wanted.")
            break
        to_outfile.append(f"Iteration {j+1}\n")
        # outfile.write(f"Iteration {j+1}\n")
        tmp_out = ""
        for i, structure in enumerate(structure_list):
            cmd.select(align_structure_list[i], structure+selection[i])
            if i != 0:
                super = cmd.super(target=align_structure_list[0], mobile=align_structure_list[i])
                print(super)
                tmp_out += f"{structure}\t{super[0]}\t{super[1]}\n"
                # outfile.write(f"{structure}\t{super[0]}\t{super[1]}\n")
                # outfile.write(f"{structure}\t{super['RMSD']}\t{super['alignment_length']}\n")
        to_outfile.append(tmp_out)
        if to_outfile[-1] == to_outfile[-3]:
            break

    with open("log.txt","w") as outfile:
        outfile.write("".join(to_outfile))

    if break_flag == False:
        print(f"Colored after {iterations} iteration(s) of superexposion.")
    for i, model in enumerate(model_kd.keys()):
        color_pymol(model, structure_list_entire[i], score_list[i])

# ### GAMMELT

#     #score_selection
#     for i, model in enumerate(models):
#         s = scores_by_conservation(model, models[:i]+models[i+1:], 6)
#         if len(s) < tresshold_aa:
#             color_pymol(model, structure_list_entire[i])
#             break_flag = True
#         elif iterations == j:
#             color_pymol(model, structure_list_entire[i], s)
#         else:
#             cmd.select(align_structure_list[i], structure_list[i]+select_by_score(s, model.atom))
#             if i != 0:
#                 try:
#                     cmd.super(align_structure_list[i], align_structure_list[0])
#                 except:
#                     cmd.super(structure_list[i], align_structure_list[0])
#     if break_flag:
#         print(f"Error! Ran {j} full iteration(s). Try setting the maximum iterations to {j}")
#         break
# if break_flag == False:
#     print(f"Ran succesfully {iterations} iteration(s)")


    #Some pymol stuff
    cmd.hide("cgo", "aln")
    cmd.set("seq_view_label_mode", "1")
    cmd.set("antialias", "4")
    cmd.set("ray_trace_mode", "1")
    # cmd.set("seq_view_fill_color", "black")

    #cmd.util.cnc("all",_self=cmd)
    print(outfile)
    cmd.save(str(outfile))

elapsed_time = time.time() - start_time
print("Elapsed time:", elapsed_time, "seconds")
# for i, model in enumerate(models): print(scores_by_conservation(model, models[:i]+models[i+1:], 6))
#super
# for i, structure in enumerate(structure_list[1:]):
#     cmd.super(structure,ref_structure)

# #udtræk models
# models = []
# for structure in structure_list:
#     models.append(cmd.get_model(structure))

#color by score



# for structure in structure_list:
#     cmd.select(structure, structure+select_by_score())
# select_string = select_by_score(s, 0.8)


# for i, model in enumerate(models):
#     s, h = scores_by_conservation(model,  models[:i]+models[i+1:], 6)
#     print(s)
# # sc, homo = scores_by_conservation(ref_model, homologous_models, 6)

# print(sc)
# su = []
# for structure in structure_list[1:]:
#     su.append(cmd.super(structure, "crystal"))

# # su = cmd.super(structure_list[1], "crystal")
# print(su)
# for i, model in enumerate(models):
#     s, h = scores_by_conservation(model,  models[:i]+models[i+1:], 6)
#     print(s)

# kage = []
# for e in structure_list:
#     kage.append(cmd.get_model(e))

# print(scores_by_conservation(kage[1], kage[2:], 6))

# # print(homo)
# print(homo[:][1][:])

# for i, model in enumerate(models): 
#     print(scores_by_conservation(model, models[:i]+models[i+1:], 6))

# make selections

# super of selection



# path_folder = "homologues"
# path_ref_structure = "ref_structure.pdb"



# ref_structure = path_ref_structure.split(".")[0]
# cmd.load(path_ref_structure)
# for file in os.listdir(path_folder): (cmd.load(path_folder + "/" + file))
# homolog_structures = cmd.get_object_list()


# structure_list = cmd.get_object_list()
# hide everything, ref_structure and chain B
# select crystal, ref_structure and resi 1:439 and chain A
# show lines, "all"




# #Blosum for alignments
# #@markdown From https://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
# def aa_to_blosom(aa1,aa2):
#     blosom62 = [
#     [4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0,-2,-1,0,-4],
#     [-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3,-1,0,-1,-4],
#     [-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3,3,0,-1,-4],
#     [-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3,4,1,-1,-4],
#     [0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4],
#     [-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2,0,3,-1,-4],
#     [-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2,1,4,-1,-4],
#     [0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3,-1,-2,-1,-4],
#     [-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3,0,0,-1,-4],
#     [-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3,-3,-3,-1,-4],
#     [-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1,-4,-3,-1,-4],
#     [-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2,0,1,-1,-4],
#     [-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1,-3,-1,-1,-4],
#     [-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1,-3,-3,-1,-4],
#     [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2,-2,-1,-2,-4],
#     [1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2,0,0,0,-4],
#     [0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0,-1,-1,0,-4],
#     [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3,-4,-3,-2,-4],
#     [-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1,-3,-2,-1,-4],
#     [0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4,-3,-2,-1,-4],
#     [-2,-1,3,4,-3,0,1,-1,0,-3,-4,0,-3,-3,-2,0,-1,-4,-3,-3,4,1,-1,-4],
#     [-1,0,0,1,-3,3,4,-2,0,-3,-3,1,-1,-3,-1,0,-1,-3,-2,-2,1,4,-1,-4],
#     [0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2,0,0,-2,-1,-1,-1,-1,-1,-4],
#     [-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,1]]
#     aa_id_list = ["ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","B","Z","X","-"]
#     aa_to_id = dict()
#     for i in range(len(aa_id_list)):
#         aa_to_id[aa_id_list[i]] = i
#     id1 = aa_to_id[aa1]
#     id2 = aa_to_id[aa2]
#     return blosom62[id1][id2]
# cmd.extend('aa_to_blosom',aa_to_blosom)
# print("Imported functions:")
# print("aa_to_blosom")

# def color_by_number(number, minmax=[0,1]):
#     """
#     Rainbow from red to blue
#     """
#     min = minmax[0]
#     max = minmax[1]
#     if number < 0.25:
#         return (max, number*4*(max-min)+min, min)
#     elif number < 0.5:
#         return ((0.25-(number-0.25))*(max-min)+min, max, min)
#     elif number < 0.75:
#         return (min, max, (number-0.5)*(max-min)+min)
#     else:
#         return (min, (0.25-(number-0.75))*(max-min)+min, max)
    
# # def color_by_number(number):
# #     """
# #     Rainbow from red to blue
# #     """
# #     if number < 0.25:
# #         return [0.9, number*1.6+0.5, 0.5]
# #     elif number < 0.5:
# #         return [1-(number-0.25)*4, 0.9, 0.5]
# #     elif number < 0.75:
# #         return [0.5, 0.9, (number-0.5)*4]
# #     else:
# #         return [0.5, 1-(number-0.75)*4, 0.9]



# # def color_by_similarity(ref_structure, homolog_list, max_dist):
# #     """Coloring ref_structure acording to sequence similarity to homolog_list with a maximum of max_dist Å between CA.
# #     blue: conserved
# #     Red: unconserved"""
# #     ref_model = cmd.get_model(f"{ref_structure} and name CA")
# #     for i, structure in enumerate(homolog_list):
# #        homolog_list[i] = cmd.get_model(f"{structure} and name CA")
# #     list_len = len(homolog_list)

# #     for ref_atom in ref_model.atom:
# #         number = 0
# #         ref_resn = ref_atom.resn
# #         ref_index = np.array(ref_atom.coord)
# #         max_number = aa_to_blosom(ref_resn,ref_resn) + 4
# #         for model in homolog_list:
# #             dist_resn = [float("inf"),""]
# #             for atom in model.atom:
# #                 atom_index = np.array(atom.coord)
# #                 dist = np.sqrt(np.sum((atom_index-ref_index)**2))
                
# #                 if dist < dist_resn[0]:
# #                     dist_resn = [dist, atom.resn]
# #             if dist_resn[0] <= max_dist:
# #                 number += ((aa_to_blosom(ref_resn,dist_resn[1]) + 4)/(list_len*max_number))
# #         if number == 0:
# #             cmd.set_color(f"color{ref_atom}", [0.5, 0.5, 0.5])
# #         elif number < 0.5:
# #             cmd.set_color(f"color{ref_atom}", [1, number*2, 0])
# #         else:
# #             cmd.set_color(f"color{ref_atom}", [1, 1, (number-0.5)*2])
# #         cmd.color(f"color{ref_atom}", f"resi {ref_atom.resi} and {ref_structure}")
# # cmd.extend('color_by_similarity',color_by_similarity)
# # print("Imported functions:")
# # print("color_by_similarity")




def make_gradient(width=10, height=100, outfile='gradient.png'):
    import png

    img = []
    for y in range(height):
        row = ()
        for x in range(width):
            number = y/100
            # if number < 0.25:
            #     row = row + (0.9, number*1.6+0.5, 0.5)
            # elif number < 0.5:
            #     row = row + (1-(number-0.25)*1.6+0.5, 0.9, 0.5)
            # elif number < 0.75:
            #     row = row + (0.5, 0.9, (number-0.5)*1.6+0.5)
            # else:
            #     row = row + (0.5, 1-(number-0.75)*1.6+0.5, 0.9)
            # kage = []
            for k in color_by_number(number):
                row += tuple([int(k*255)])
                # kage.append(int(k*255))
            # row = row + (kage)
        img.append(row)
    with open(outfile, 'wb') as f:
        w = png.Writer(width, height, greyscale=False)
        w.write(f, img)
    for x in range(10):
        print(color_by_number(x/10)*255)

# import png

# width = 255
# height = 255
# img = []
# for y in range(height):
#     row = ()
#     for x in range(width):
#         row = row + (x, max(0, 255 - x - y), y)
#     img.append(row)
# print(img)
# with open('gradient.png', 'wb') as f:
#     w = png.Writer(width, height, greyscale=False)
#     w.write(f, img)



# # Main code - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -




# # make selections

# # super of selection


# # for structure in homolog_structures: cmd.super(structure, "crystal")
# for structure in homolog_structures: cmd.cealign(target="crystal", mobile=structure, object="aln")
# save alignment.aln, aln
# for structure in homolog_structures: cmd.super(structure, "crystal")

# run color_by_similarity_v2.py

# #color all structures by similarity
# for i, structure in enumerate(structure_list): color_by_similarity(structure, structure_list[:i]+structure_list[i+1:], 6)

# #color only SP3 by similarity
# #color_by_similarity("crystal", homolog_structures, 6)
# util.cnc("all",_self=cmd)





# # ref_structure = path_ref_structure.split(".")[0]
# # cmd.load(path_ref_structure)
# # for file in os.listdir(path_folder): (cmd.load(path_folder + "/" + file))
# # homolog_structures = cmd.get_object_list()


# # structure_list = cmd.get_object_list()
# # hide everything, ref_structure and chain B
# # select crystal, ref_structure and resi 1:439 and chain A
# # show lines, "all"