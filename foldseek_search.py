# def get_file(file):
#     raw_script = f"https://raw.githubusercontent.com/RuneROe/git_color_by_similarity/master/{file}"
#     local_script_path = f"/content/{script}.py"
#     os.system(f"wget {raw_script} -O {local_script_path}")


def run(database,variable_tresshold,value_tresshold,search_against,ref_structure,infilenames): 
    low_is_good = {"evalue","rmsd"}   
    if database == "Thermophilic_DB":
        print("Downloading thermophilic database")
        ##################Fix
        DB = "thermoDB/thermoDB"
    else:
        print("Downloading database:",database)
        os.system(f"foldseek databases {database} DB tmp")
        os.system(f"foldseek createindex DB tmp")
        DB = "DB"
    os.system("mkdir foldseek_output")
    if search_against == "ref_structure":
        os.system(f"foldseek easy-search {ref_structure} {DB} foldseek_output/aln.txt tmp --format-output target,{variable_tresshold}")
    else:
        for structure in infilenames:
            os.system(f"foldseek easy-search {structure} {DB} foldseek_output/{structure.split(".")[0]}.txt tmp --format-output target,{variable_tresshold}")