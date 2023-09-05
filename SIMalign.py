
#TO DO
#- Lav bedre output fra run function
#- rydde op i py doc
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

    for i, atom in enumerate(model_CA.atom):
        cmd.set_color(f"{atom}color", color_by_number(score_list[i]))
        cmd.color(f"{atom}color", f"resi {atom.resi} and chain {atom.chain} and {structure_entire}")



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
    




def run(ref_structure, files, iterations, tresshold_aa, max_dist, remove_chain_duplicate, outfilename):
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


# Importing files - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - -
    try:
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
        structure_list = []
        align_structure_list = []
        for i, structure in enumerate(structure_list_entire):
            structure_list.append(structure+"_CA")
            align_structure_list.append(structure+"_align")
            cmd.select(structure_list[i], structure+" and name CA and not HETATM")
    except:
        print("\n- - - - - - - - - - - - - - - - - -\nImport ERROR")
        print("\nCouldn't import one of the files in pymol\n- - - - - - - - - - - - - - - - - -\n")
        sys.exit(1)

# Basic Pymol stuff - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - -

    cmd.remove("hydrogens")
    cmd.alignto(ref_structure, object="aln")

# LOOP start
    break_flag = False
    to_outfile = []
    for j in range(iterations):
        model_kd = dict()  
        for structure in structure_list:
            model = cmd.get_model(structure)
            model_kd[model] = cKDTree([atom.coord for atom in model.atom])
        score_list = []
        selection = []
        for i, model in enumerate(model_kd.keys()):
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
                        if dist_resn[1] <= max_dist:
                            s += ((aa_to_blosom(ref_resn,dist_resn[2]) + 4)/(n_homologous_list*max_score))
                            tmp.append(dist_resn[1:])
                        dist_list.append(dist_resn)
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
        to_outfile.append(f"Iteration {j+1}\nstructure\tRMSD\tatoms aligned\n")
        tmp_out = ""
        for i, structure in enumerate(structure_list):
            cmd.select(align_structure_list[i], structure+selection[i])
            if i != 0:
                super = cmd.super(target=align_structure_list[0], mobile=align_structure_list[i])
                tmp_out += f"{structure_list_entire[i]}\t{round(super[0],3)}\t{super[1]}\n"
        to_outfile.append(tmp_out)
        print(to_outfile[-2]+to_outfile[-1])
        if to_outfile[-1] == to_outfile[-3]:
            break
        
    with open("log.txt","w") as outfile:
        outfile.write("".join(to_outfile))

    if break_flag == False:
        print(f"Colored after {iterations} iteration(s) of superexposion.")
    for i, model in enumerate(model_kd.keys()):
        color_pymol(model, structure_list_entire[i], score_list[i])


# Some pymol stuff - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - -
    cmd.hide("cgo", "aln")
    cmd.set("seq_view_label_mode", "1")
    cmd.set("antialias", "4")
    cmd.set("ray_trace_mode", "1")

    cmd.save(outfilename)



def make_gradient(width=10, height=100, outfile='gradient.png'):
    import png

    img = []
    for y in range(height):
        row = ()
        for x in range(width):
            number = y/100
            for k in color_by_number(number):
                row += tuple([int(k*255)])
        img.append(row)
    with open(outfile, 'wb') as f:
        w = png.Writer(width, height, greyscale=False)
        w.write(f, img)

