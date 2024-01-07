import os
mainfolder = "homologs_foldseekAF-proteome"
for file in os.listdir(mainfolder):
	file = f"{mainfolder}/{file}"
	folder = file.split(".")[0]
	print(f"Making folder: {folder}")
	os.system(f"mkdir {folder}")
	with open(file,"r") as infile:
		lines = infile.readlines()[:20]
		for line in lines:
			if line.endswith("cif.gz"):
				cif = ".".join(line.split("\t")[1].split(".")[:-1])
			else:
				cif = line.split("\t")[1]+".cif"
			print(f"Downloading: {cif} to {folder}")
			os.system(f"gsutil -m cp gs://public-datasets-deepmind-alphafold-v4/{cif} {folder}/.")