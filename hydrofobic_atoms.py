from pymol import cmd
from scipy.spatial import cKDTree

print("Imported functions:")
cmd.reinitialize()

def color_by_number(number):
    """
    DESCRIPTION

    Rainbow from red to white
    red white
    number 0: blue
    number 0.5: white
    number 1: red

    DEPENDENCIES

    import numpy as np
    """
    if number < 0.5:
        return [0.2+(1-0.2)*number*2,0.937+(1-0.937)*number*2,0.44+(1-0.44)*number*2]
    else:
        return [0.937+(1-0.937)*(1-(number-0.5)*2),0.44+(1-0.44)*(1-(number-0.5)*2),0.2+(1-0.2)*(1-(number-0.5)*2)]
    if number >= 1:
        return [1,1,1]
    return [0.8+(number/5),number,number]





def name_to_hydrofobic(atom_name):
    logD_by_aa = {# https://web.expasy.org/protscale/pscale/Hphob.Eisenberg.html
    "ALA":  0.620,  
    "ARG": -2.530,  
    "ASN": -0.780,  
    "ASP": -0.900,  
    "CYS":  0.290,  
    "GLN": -0.850,  
    "GLU": -0.740,  
    "GLY":  0.480,  
    "HIS": -0.400,  
    "ILE":  1.380,  
    "LEU":  1.060,  
    "LYS": -1.500,  
    "MET":  0.640,  
    "PHE":  1.190,  
    "PRO":  0.120,  
    "SER": -0.180,  
    "THR": -0.050,  
    "TRP":  0.810,  
    "TYR":  0.260,  
    "VAL":  1.080}
    new = {}
    for k,v in logD_by_aa.items():
        new[k] = v
    side = {
    "ALA":  1,  
    "ARG": 10,  
    "ASN": 4,  
    "ASP": 4,  
    "CYS":  2,  
    "GLN": 5,  
    "GLU": 5,  
    "GLY":  0,  
    "HIS": 7,  
    "ILE":  4,  
    "LEU":  4,  
    "LYS": 5,  
    "MET":  4,  
    "PHE":  6,  
    "PRO": 3,  
    "SER": 2,  
    "THR": 3,  
    "TRP":  10,  
    "TYR":  7,  
    "VAL":  3}
    return new[atom_name]


def hydrofobic_side_chain(selection="all"):
    model = cmd.get_model(selection+" and not HETATM").atom
    # side_chain = cmd.get_model(selection+" and (not name O and not name CA and not name C and not name N)")
    KD = cKDTree([atom.coord for atom in model])
    score = {}
    for atom in model:
        near = KD.query_ball_point(atom.coord,4)
        hydro_value = 0
        len_near = 0
        for x in near:
            if model[x].name not in {"CA","C","N","O"}:
                hydro_value += name_to_hydrofobic(model[x].resn)
            len_near += 1
        if len_near != 0:
            hydro_value = hydro_value/len_near
        else:
            print(selection,atom.resi,atom.name,[model[x].name for x in near])
        # if hydro_value < 0:
        #     hydro_value = (hydro_value + 3.01)/(2*3.01) 
        # else:
        #     hydro_value = ((hydro_value/0.9)/2)+0.5 
        if hydro_value < 0:
            hydro_value = (hydro_value + 2.53)/(2*2.53) #2.53
        else:
            hydro_value = ((hydro_value/1.38)/2)+0.5 #1.38
        # hydro_value = (hydro_value+2.53)/3.91
        # hydro_value = 0
        # if atom.name not in {"CA","C","N","O"}:
        #     hydro_value = name_to_hydrofobic(atom.resn)
        # hydro_value = (hydro_value+3.01)/3.91
        cmd.set_color(atom,color_by_number(hydro_value))
        cmd.color(atom,selection+f" and chain {atom.chain} and resi {atom.resi} and name {atom.name}")
        score[atom.resi] = hydro_value
    return score
cmd.extend('hydrofobic_side_chain',hydrofobic_side_chain)
print("hydrofobic_side_chain")





cmd.load("../SP/SP3.pdb")
# cmd.load("../SP/SP3_apo.pdb")
# cmd.load("../SP/SP1.pdb")
# cmd.load("../SP/SP2.pdb")
cmd.hide("everything","chain B")
# cmd.alignto("SP3")
structures = cmd.get_object_list()
cmd.remove("hydrogens")
for s in structures:
    x = hydrofobic_side_chain(s+" and chain A")
cmd.show("spheres","chain A")
cmd.save("tmp/h_test.pse")
print(x)


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

make_gradient(outfile="hydro_gradient.png")