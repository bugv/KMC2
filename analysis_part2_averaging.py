import numpy as np
import h5py 
import os 
import copy

os.makedirs("results_second_processing", exist_ok=True)

species = ["Al","Fe","X0+"]

for root,subdir,files in os.walk("results_first_processing"):
    if len(files) > 100 :
        print("Processing", root)
        print("hi")
        os.makedirs(os.path.join("results_second_processing",*root.split("/")[1:-1]), exist_ok=True)
        output_filename = os.path.join("results_second_processing",*root.split("/")[1:-1], root.split("/")[-1] + ".h5")

        with h5py.File(output_filename,"w") as outfile :

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


            for key in ["r2s", "theproduct"] :
                data_group_mean = h5file.create_group(key + "_mean")
                data_group_std = h5file.create_group(key + "_std")
                for subkey, subvalue in value.items():
                    data_group_mean.create_dataset(subkey, data=subvalue.mean(axis=1))
                    data_group_std.create_dataset(subkey, data=subvalue.std(axis=1))

            outfile.create_dataset("time_collector", data=time_collector.mean(axis=1))

            outfile.create_dataset("time_collector_std", data=time_collector.std(axis=1))

            outfile.create_dataset("n_runs", data=n_runs)