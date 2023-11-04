import py3Dmol
import SIMalign

view = py3Dmol.view(js='https://3dmol.org/build/3Dmol.js',)
view.addModel(open(ref_structure,'r').read(),'pdb')

if color_mode == "similarity":
    color_list = []
    for score in score_list[0]:
        tmp = SIMalign.color_by_number(score)
        string = "rgb("
        for number in tmp:
            string += str(int(number*255))+","
        string[:-1]+")"
        color_list.append(string)

    view.setStyle({'cartoon':{'colorscheme':{'prop':'resi',"gradient":"linear",'colors':color_list,'min':1,'max':len_ref_structure}}})
elif color_mode == "hotspot":
    pass
else:
    view.setStyle({'cartoon': {'color':'spectrum'}})


view.zoomTo()