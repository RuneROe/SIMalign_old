def color_by_number(number):
    """
    DESCRIPTION

    Rainbow from red to white
    number 0: red
    number 1: white

    DEPENDENCIES

    import numpy as np
    """
    if number >= 1:
        return [1,1,1]
    return [0.8+(number/5),number,number]

from pymol import cmd

# color_by_score(structure, score_list[j])




def list_to_selection(list,select_above99=False):
    selection_list = []
    for ele in list:
        selection = " and ("
        if isinstance(ele, dict):
            ele = ele.keys()
        if select_above99:
            for i, score in enumerate(ele):
                if float(score)>0.99:
                    selection += f"resi {i+1} or "
        else:
            for resi in ele:
                selection +=  f"resi {int(resi)} or "
        if selection == " and (":
            selection_list.append("")
        else:
            selection_list.append(selection[:-4]+")")
    return selection_list




def select_by_list(selection_list,structure_list,list_of_lists=True,select_above99=False):
    if list_of_lists and select_above99:
        selection_list = list_to_selection(selection_list,select_above99=True)
    elif list_of_lists:
        selection_list = list_to_selection(selection_list)
    selection = " and (("
    for i, structure in enumerate(structure_list):
        if selection_list[i] != "":
            selection += structure+selection_list[i]+") or ("
    if len(selection) < 10:
        return ""
    else:
        return selection[:-4]+")"


def color_by_score(structure, score):
    """
    DESCRIPTION



    DEPENDENCIES

    color_by_number
    from pymol import cmd
    """
    model = cmd.get_model(structure)
    resi1 = int(model.atom[0].resi)
    # select_conserved_list = []
    for i, s in enumerate(score):
        cmd.set_color(f"{str(model.atom[i])+structure}color", color_by_number(s))
        cmd.color(f"{str(model.atom[i])+structure}color", f"resi {resi1+i} and {structure}")
        # if score[i] > 0.99:
            # select_conserved_list.append(int(atom.resi))
    # cmd.select(f"{structure}_conserved",structure+f" and not HETATM and chain A{list_to_selection(select_conserved_list)}")




# def color_by_hotspot(structure, hotspot):
#     grey = [0.9,0.9,0.9]
#     colors = [
#         [0.81, 0.34, 0.34], #red
#         [0.42, 0.47, 0.91], #blue
#         [0.4, 0.84, 0.47], #lime
#         [0.82, 0.38, 0.83], #purble
#         [0.91, 0.64, 0.45], #
#         [0.41, 0.89, 0.7], 
#         [0.27, 0.71, 0.84], 
#         [0.69, 0.56, 0.82], 
#         [0.9, 0.68, 0.77], 
#         [0.88, 0.79, 0.42]]
#     model = cmd.get_model(structure+" and name CA and not HETATM and chain A")
#     for i, atom in enumerate(model.atom):
#         not_in_hotspot = True
#         for j, h in enumerate(hotspot):
#             cmd.select(f"{structure}_HS{j+1}",f"{structure} and not HETATM and chain A{list_to_selection(h)}")
#             if i in h:
#                 not_in_hotspot = False
#                 cmd.set_color(f"{atom}color", colors[j%10])
#                 cmd.color(f"{atom}color", f"resi {atom.resi} and chain {atom.chain} and {structure}")
#         if not_in_hotspot:
#             cmd.set_color(f"{atom}color", grey)
#             cmd.color(f"{atom}color", f"resi {atom.resi} and chain {atom.chain} and {structure}")


# def simple_color_by_hotspot(structure, hotspot):
#     cmd.set_color("hot",[0.81, 0.34, 0.34])
#     cmd.set_color("grey",[0.9,0.9,0.9])
#     cmd.color("grey", structure+" and chain A and not HETATM")
#     if len(hotspot) > 0:
#         cmd.select(f"{structure}_HS",f"{structure} and not HETATM and chain A{list_to_selection(hotspot.keys())}")
#         cmd.color("hot", structure+"_HS")



def structure_formating(structure_format,show_in_pymol):
    cmd.hide("everything","not chain A")
    if structure_format == "cartoon":
        cmd.show_as("cartoon","not HETATM")
    elif structure_format == "spheres":
        cmd.show_as("spheres","not HETATM")
    else:
        cmd.show_as("sticks","not HETATM")
        if structure_format == "spheres-sticks":
            cmd.show("spheres","not HETATM")
            cmd.set("sphere_transparency", "0.7")
    if show_in_pymol != "everything":
        cmd.hide("everything","not chain A")
        if show_in_pymol != "entire_chain_A" and show_in_pymol != "only_not_conserved":
            cmd.hide("everything","exposed_AA")
            if show_in_pymol != "only_not_conserved_core":
                cmd.hide("everthing","conserved")
        elif show_in_pymol == "only_not_conserved":
            cmd.hide("everything","conserved")


def select_exposed_AA(structure_list):
    import findsurfaceatoms
    exposed_list = []
    for structure in structure_list:
        exposed_residues = findsurfaceatoms.findSurfaceResidues(selection=structure)
        exposed_set = set()
        for resi in exposed_residues:
            exposed_set.add(resi[1])
        exposed_list.append(exposed_set)
    cmd.select("exposed_AA","not HETATM"+select_by_list(exposed_list,structure_list))

import SIMalign

def run(color_mode,hotspot_list,score_list,structure_list,core_selection,exposed_list,structure_format,show_in_pymol,color_by_element):
    print("Coloring and selecting in pymol...")
    structure_list_entire_chainA = SIMalign.select_first_chain(structure_list)
    #Select core_AA
    # print("not HETATM"+select_by_list(core_selection,structure_list,list_of_lists=False))
    # print("not HETATM and chain A"+select_by_list(score_list,structure_list,select_above99=True))
    # print("not conserved and not HETATM and chain A")
    # print("chain A and not HETATM"+select_by_list(hotspot_list,structure_list))
    # print("chain A and not HETATM"+select_by_list(exposed_list,structure_list))
    cmd.select("core_AA","not HETATM"+select_by_list(core_selection,structure_list_entire_chainA,list_of_lists=False))

    #Select conserved and nonconserved
    cmd.select("conserved", "not HETATM"+select_by_list(score_list,structure_list_entire_chainA,select_above99=True))
    cmd.select("nonconserved","not conserved and not HETATM and chain A")
    
    if hotspot_list != None:
        #Select hotspots and exposed AA
        cmd.select("hotspots","not HETATM"+select_by_list(hotspot_list,structure_list_entire_chainA))
    if exposed_list != None and exposed_list != []:
        cmd.select("exposed_AA","not HETATM"+select_by_list(exposed_list,structure_list_entire_chainA))
    else:
        select_exposed_AA(structure_list_entire_chainA)
    structure_formating(structure_format,show_in_pymol)
    if color_mode == "hotspot":
        print("\tColoring by hotspots")
        cmd.set_color("hot",[0.82, 0.38, 0.83])
        cmd.set_color("grey",[0.9,0.9,0.9])
        cmd.color("grey","not HETATM and not hotspots")
        cmd.color("hot","hotspots")
        # for j, structure in enumerate(structure_list):
        #     print(f"Coloring {structure}")
        #     simple_color_by_hotspot(structure, hotspot_list[j])
    elif color_mode == "similarity":
        print("\tColoring by similarity:")
        for j, structure in enumerate(structure_list_entire_chainA):
            print(f"\t\tColoring {structure_list_entire_chainA[j].split(' ')[0]}")
            color_by_score(structure, score_list[j])
    if color_by_element:
        cmd.util.cnc("all")
    cmd.hide("cgo", "aln")
    cmd.set("seq_view_label_mode", "1")
    cmd.set("antialias", "4")
    cmd.set("ray_trace_mode", "1")
    # cmd.save(outfile_name)

# colors.run("similarity",hotspot_list,score_list,outfile_name,structure_list)

# def make_gradient(width=10, height=100, outfile='gradient.png'):
#     import png

#     img = []
#     for y in range(height):
#         row = ()
#         for x in range(width):
#             number = y/100
#             for k in color_by_number(number):
#                 row += tuple([int(k*255)])
#         img.append(row)
#     with open(outfile, 'wb') as f:
#         w = png.Writer(width, height, greyscale=False)
#         w.write(f, img)


import py3Dmol

def show_pdb(ref_structure,color_mode,score_list,len_ref_structure,hotspot_list):
    view = py3Dmol.view(js='https://3dmol.org/build/3Dmol.js',)
    view.addModel(open(ref_structure,'r').read(),'pdb')

    if color_mode == "similarity":
        color_list = []
        for score in score_list[0]:
            tmp = color_by_number(score)
            string = "rgb("
            for number in tmp:
                string += str(int(number*255))+","
            string[:-1]+")"
            color_list.append(string)

        view.setStyle({'cartoon':{'colorscheme':{'prop':'resi',"gradient":"linear",'colors':color_list,'min':1,'max':len_ref_structure}}})
    elif color_mode == "hotspot":
        color_list = []
        h = set(hotspot_list[0].keys())
        for index in range(len(score_list[0])):
            string = "rgb()"
            if index+1 in h:
                string = "rgb(85,225,175)"
            else:
                string = "rgb(220,220,220)"
            color_list.append(string)
        view.setStyle({'cartoon':{'colorscheme':{'prop':'resi',"gradient":"linear",'colors':color_list,'min':1,'max':len_ref_structure}}})
    else:
        view.setStyle({'cartoon': {'color':'spectrum'}})
    return view

    