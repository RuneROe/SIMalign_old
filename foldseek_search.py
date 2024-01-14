import os
import urllib.request

def download_AF(name,outfolder,infilenames):
    if name.endswith("cif.gz"):
        pdb = ".".join(name.split(".")[:-2])+".pdb"
    else:
        pdb = name+".pdb"
    print("\tDownloading:",pdb)
    url = "https://alphafold.ebi.ac.uk/files/"+pdb
    urllib.request.urlretrieve(url,outfolder+"/"+pdb)
    infilenames.append(outfolder+"/"+pdb)
    return infilenames

def foldseek_search(structure, DB,variable_tresshold):
    if variable_tresshold == "number_of_structures":
        os.system(f"./foldseek/bin/foldseek easy-search {structure} {DB} foldseek_output/{structure.split('.')[0]}.txt tmp --format-output target")
    else:
        os.system(f"./foldseek/bin/foldseek easy-search {structure} {DB} foldseek_output/{structure.split('.')[0]}.txt tmp --format-output target,{variable_tresshold}")


def remove_old_structures():
    if "foldseek_output" in os.listdir():
        if "structures" in os.listdir("foldseek_output"):
            for file in os.listdir("foldseek_output/structures"):
                os.remove("foldseek_output/structures/"+file)

def run(database,variable_tresshold,value_tresshold,search_against,ref_structure,infilenames): 
    print("Running foldseek...")
    low_flag = False
    low_is_good = {"evalue","rmsd"}   
    if variable_tresshold in low_is_good:
        low_flag = True
    print("mayby thermo")
    if database == "Thermophilic_DB":
        print("thermo")
        if not os.path.isfile("ThermoDB_READY"): 
            print("\tDownloading thermophilic database")
            os.system("pip install gdown")
            print("\t\tgdown installed")
            # os.system("gdown --folder https://drive.google.com/drive/u/1/folders/1FN3Cfl94J0ML2UmRADNFuTAqabOkxdfN")
            os.system("gdown https://drive.google.com/file/d/1255hcwjyDTE7tR7Ahmik57W7cgyYuHkP/view?usp=sharing")
            print("\t\tThermophilic DB downloaded")
            os.system("touch ThermoDB_READY")
        DB = "thermoDB/thermoDB"
    else:
        DB = "DB"+database.split("/")[-1]
        if not os.path.isfile("foldseek_"+DB):
            print("\tDownloading database:",database)
            os.system(f"./foldseek/bin/foldseek databases {database} {DB} tmp")
            os.system(f"./foldseek/bin/foldseek createindex {DB} tmp")
            os.system("touch foldseek_"+DB)
    remove_old_structures()
    os.system("mkdir foldseek_output")
    if search_against == "ref_structure":
        foldseek_search(ref_structure, DB,variable_tresshold)
    else:
        for structure in infilenames:
            foldseek_search(structure, DB,variable_tresshold)
    if "structures" in os.listdir("foldseek_output"):
        for file in os.listdir("foldseek_output/structures"):
            os.remove("foldseek_output/structures/"+file)
    else:
        os.system("mkdir foldseek_output/structures")
    for file in os.listdir("foldseek_output"):
        if file.endswith(".txt"):
            with open("foldseek_output/"+file) as infile:
                if variable_tresshold == "number_of_structures":
                    lines = infile.readlines()
                    for i in range(int(value_tresshold)):
                        infilenames = download_AF(lines[i][:-1],"foldseek_output/structures",infilenames)
                else:
                    for line in infile:
                        line_list = line.split("\t")
                        variable = float(line_list[1][:-1])
                        if low_flag:
                            if variable < value_tresshold:
                                infilenames = download_AF(line_list[0],"foldseek_output/structures",infilenames)
                        elif variable > value_tresshold:
                            infilenames = download_AF(line_list[0],"foldseek_output/structures",infilenames)
    return infilenames