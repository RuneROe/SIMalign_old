# def get_file(file):
#     raw_script = f"https://raw.githubusercontent.com/RuneROe/git_color_by_similarity/master/{file}"
#     local_script_path = f"/content/{script}.py"
#     os.system(f"wget {raw_script} -O {local_script_path}")
import os

def download_AF(name,outfolder,infilenames):
    if name.endswith("cif.gz"):
        cif = name.split(".")[:-1]
    else:
        cif = name+".cif"
    print("/tDownloading:",cif)
    os.system(f"gsutil -m cp gs://public-datasets-deepmind-alphafold-v4/{cif} {outfolder}/.")
    infilenames.append(outfolder+"/"+cif)
    return infilenames


def run(database,variable_tresshold,value_tresshold,search_against,ref_structure,infilenames): 
    print("Running foldseek...")
    low_flag = False
    low_is_good = {"evalue","rmsd"}   
    if variable_tresshold in low_is_good:
        low_flag = True
    if database == "Thermophilic_DB":
        if not os.path.isfile("ThermoDB_READY"): 
            print("Downloading thermophilic database")
            os.system("pip install gdown")
            os.system("gdown --folder https://drive.google.com/drive/u/1/folders/1FN3Cfl94J0ML2UmRADNFuTAqabOkxdfN")
            os.system("touch ThermoDB_READY")
        DB = "thermoDB/thermoDB"
    else:
        DB = "DB"+database.split("/")[-1]
        if not os.path.isfile("foldseek_"+DB):
            print("Downloading database:",database)
            os.system(f"./foldseek/bin/foldseek databases {database} {DB} tmp")
            os.system(f"./foldseek/bin/foldseek createindex {DB} tmp")
            os.system("touch foldseek_"+DB)
    os.system("mkdir foldseek_output")
    if search_against == "ref_structure":
        os.system(f"./foldseek/bin/foldseek easy-search {ref_structure} {DB} foldseek_output/aln.txt tmp --format-output target,{variable_tresshold}")
    else:
        for structure in infilenames:
            os.system(f"./foldseek/bin/foldseek easy-search {structure} {DB} foldseek_output/{structure.split('.')[0]}.txt tmp --format-output target,{variable_tresshold}")
    os.system("mkdir foldseek_output/structures")
    for file in os.listdir("foldseek_output"):
        if file != "structures":
            with open("foldseek_output/"+file) as infile:
                for line in infile:
                    line_list = line.split("/t")
                    variable = line_list[1]
                    if low_flag:
                        if variable < value_tresshold:
                            infilenames = download_AF(line_list[0],"foldseek_output/structures",infilenames)
                    elif variable > value_tresshold:
                        infilenames = download_AF(line_list[0],"foldseek_output/structures",infilenames)
    return infilenames


    
    

# import os
# mainfolder = "homologs_foldseekAF-proteome"
# for file in os.listdir(mainfolder):
# 	file = f"{mainfolder}/{file}"
# 	folder = file.split(".")[0]
# 	print(f"Making folder: {folder}")
# 	os.system(f"mkdir {folder}")
# 	with open(file,"r") as infile:
# 		lines = infile.readlines()[:20]
# 		for line in lines:
# 			if line.endswith("cif.gz"):
# 				cif = ".".join(line.split("\t")[1].split(".")[:-1])
# 			else:
# 				cif = line.split("\t")[1]+".cif"
# 			print(f"Downloading: {cif} to {folder}")
# 			os.system(f"gsutil -m cp gs://public-datasets-deepmind-alphafold-v4/{cif} {folder}/.")