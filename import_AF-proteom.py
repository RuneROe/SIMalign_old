import os

#dependencies: gsutil

down = set()
for file in os.listdir("imported_tars"):
    if file.endswith(".tar"):
        file = file.split("tax_id-")[-1].split("-")[0]
        down.add(file)

with open("ThermoBase.tax","r") as file:
    for line in file:
        tax = line[:-1]
        if tax not in down:
            print(f"{tax}")
            try:
                os.system(f"gsutil -m cp gs://public-datasets-deepmind-alphafold-v4/proteomes/proteome-tax_id-{tax}-*_v4.tar imported_tars/.")
                down.add(tax)
            except:
                print(f"ERROR! Could not download {tax}")


with open("down_tax.txt","w") as d:
    for id in down:
        d.write(f"{id}\n")

