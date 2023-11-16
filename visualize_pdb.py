import py3Dmol
import SIMalign

def run(ref_structure,color_mode,score_list,len_ref_structure,hotspot_list):
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


    view.zoomTo()