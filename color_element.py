#     colors = [
#         [0.81, 0.34, 0.34], #red
#         [0.42, 0.47, 0.91], #blue
#         [0.4, 0.84, 0.47], #lime
#         [0.82, 0.38, 0.83], #purble
#         [0.91, 0.64, 0.45], #xbeige
#         [0.41, 0.89, 0.7], #light green
#         [0.27, 0.71, 0.84], #light blue
#         [0.69, 0.56, 0.82], #purble
#         [0.9, 0.68, 0.77], #lyserød
#         [0.88, 0.79, 0.42]]

from pymol import cmd

def color_structure(color,selection="all"):
    # cmd.set_color("red","[0.81, 0.34, 0.34]")
    # cmd.set_color("blue","[0.42, 0.47, 0.91]")
    # cmd.set_color("yellow","[0.88, 0.79, 0.42]")
    # cmd.set_color("lime","[0.4, 0.84, 0.47]")
    # cmd.set_color("purble","[0.82, 0.38, 0.83]")
    # cmd.set_color("x","[0.91, 0.64, 0.45]")
    # cmd.set_color("y","[0.41, 0.89, 0.7]")
    # cmd.set_color("z","[0.27, 0.71, 0.84]")
    # cmd.set_color("w","[0.69, 0.56, 0.82]")
    # cmd.set_color("q","[0.9, 0.68, 0.77]")
    cmd.set_color("water","[0.863, 0.365, 0.871]")
    cmd.set_color("green","[0.376, 0.827, 0.447]")
    cmd.set_color("orange","[0.933, 0.6, 0.176]")
    cmd.set_color("yellow","[0.953, 0.835, 0.349]")
    cmd.set_color("blue","[0.365, 0.416, 0.91]")
    cmd.set_color("red","[0.847, 0.322, 0.322]")
    cmd.color(color,selection)
    model = cmd.get_model(selection).atom
    for atom in model:
        if atom.name.startswith("N"):
            cmd.color("blue",selection+" and resi "+atom.resi+" and name "+atom.name)
        elif atom.name.startswith("O"):
            cmd.color("red",selection+" and resi "+atom.resi+" and name "+atom.name)
        elif atom.name.startswith("S"):
            cmd.color("yellow",selection+" and resi "+atom.resi+" and name "+atom.name)


# [0.863, 0.365, 0.871]
# [0.376, 0.827, 0.447]
# [0.933, 0.6, 0.176]
# [0.953, 0.835, 0.349]
# [0.365, 0.416, 0.91]
# [0.847, 0.322, 0.322]