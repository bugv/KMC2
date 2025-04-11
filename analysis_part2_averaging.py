import numpy as np
import h5py 
import os 
import copy
import argparse

parser = argparse.ArgumentParser(
    description="Averages of series of KMC runs, outputs of which (after processing) are given in a single directory."
)
parser.add_argument("head", type=str, help="where to start looking from")


os.makedirs("results_second_processing", exist_ok=True)

for root,subdir,files in os.walk("results_first_processing"):
    if len(files) > 100 :
        print("Processing", root)
        os.makedirs(os.path.join("results_second_processing",root.split("/")[-2]), exist_ok=True)
        output_filename = os.path.join("results_second_processing",root.split("/")[-2], root.split("/")[-1] + ".h5")

        with h5py.File(output_filename,"w") as outfile :

            if files[0].endswith(".h5") : 
                with h5py.File(os.path.join(root,files[0]),"r") as f :
                    species = list(f["input"]["atom_key"].keys())

            else : 
                    print(f"{f} is not an h5 file. Composition needs to be manually specified")

            n_runs = sum([file.endswith(".h5") for file in os.listdir(root)])

            r2s = {}
            theprod = {}
            for atom1 in species : 
                r2s[atom1] = np.full((1000,n_runs),np.nan)
                for atom2 in species :
                    if atom1 <= atom2 :
                        theprod[atom1+atom2] = np.full((1000,n_runs),np.nan)
            
            time_collector = np.full((1000,n_runs),np.nan)


            i = 0
            for file in files:
                if file.endswith(".h5"):
                    input_filename = os.path.join(root,file)
                    with h5py.File(input_filename, "r") as f:

                        for atom1 in species :
                            r2s[atom1][:, i] = f["r2s"][atom1][:]
                            for atom2 in species :
                                if atom1 <= atom2 :
                                    theprod[atom1+atom2][:, i] = f["theproduct"][atom1+atom2][:]
            
                        
                        time_collector[:, i] = f["time_collector"][:]

                        if i == 0 :
                            f.copy("input", outfile)

                        i += 1

            data_group_r2s_mean = outfile.create_group("r2s_mean")
            data_group_r2s_std = outfile.create_group("r2s_std")
            data_group_theproduct_mean = outfile.create_group("theproduct_mean")
            data_group_theproduct_std = outfile.create_group("theproduct_std")

                
            for key, value in r2s.items():
                data_group_r2s_mean.create_dataset(key, data=value.mean(axis=1))
                data_group_r2s_std.create_dataset(key, data=value.std(axis=1))
            
            for key, value in theprod.items():
                data_group_theproduct_mean.create_dataset(key, data=value.mean(axis=1))
                data_group_theproduct_std.create_dataset(key, data=value.std(axis=1))

            outfile.create_dataset("time_collector", data=time_collector.mean(axis=1))

            outfile.create_dataset("time_collector_std", data=time_collector.std(axis=1))

            outfile.create_dataset("n_runs", data=n_runs)