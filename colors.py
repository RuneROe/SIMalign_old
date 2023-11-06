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

from pymol import cmd

# color_by_score(structure, score_list[j])

def color_by_score(structure, score):
    """
    DESCRIPTION



    DEPENDENCIES

    color_by_number
    from pymol import cmd
    """
    model = cmd.get_model(structure+" and name CA and not HETATM and chain A")
    for i, atom in enumerate(model.atom):
        cmd.set_color(f"{str(atom)+structure}color", color_by_number(score[i]))
        cmd.color(f"{str(atom)+structure}color", f"resi {atom.resi} and chain {atom.chain} and {structure}")
    # cmd.save(outfile_name)

def hotspot_to_selection(hotspot):
    selection = " and ("
    for index in hotspot:
        selection +=  f"resi {index+1} or "
    return selection[:-4]+")"


def color_by_hotspot(structure, hotspot):
    grey = [0.9,0.9,0.9]
    colors = [
        [0.81, 0.34, 0.34], 
        [0.42, 0.47, 0.91], 
        [0.4, 0.84, 0.47], 
        [0.82, 0.38, 0.83], 
        [0.91, 0.64, 0.45], 
        [0.41, 0.89, 0.7], 
        [0.27, 0.71, 0.84], 
        [0.69, 0.56, 0.82], 
        [0.9, 0.68, 0.77], 
        [0.88, 0.79, 0.42]]
    model = cmd.get_model(structure+" and name CA and not HETATM and chain A")
    for i, atom in enumerate(model.atom):
        not_in_hotspot = True
        for j, h in enumerate(hotspot):
            cmd.select(f"{structure}_HS{j+1}",f"{structure} and not HETATM and chain A{hotspot_to_selection(h)}")
            if i in h:
                not_in_hotspot = False
                cmd.set_color(f"{atom}color", colors[j%10])
                cmd.color(f"{atom}color", f"resi {atom.resi} and chain {atom.chain} and {structure}")
        if not_in_hotspot:
            cmd.set_color(f"{atom}color", grey)
            cmd.color(f"{atom}color", f"resi {atom.resi} and chain {atom.chain} and {structure}")

def run(color_mode,hotspot_list,score_list,outfile_name,structure_list):
    if color_mode == "hotspot":
        print("Coloring by hotspots:\n")
        for j, structure in enumerate(structure_list):
            print(f"Coloring {structure}")
            color_by_hotspot(structure, hotspot_list[j])
    elif color_mode == "similarity":
        print("Coloring by similarity:\n")
        for j, structure in enumerate(structure_list):
            print(f"Coloring {structure}")
            color_by_score(structure, score_list[j])
    # cmd.save(outfile_name)

# colors.run("similarity",hotspot_list,score_list,outfile_name,structure_list)

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