import numpy as np
import h5py 
import os 
import copy

os.makedirs("results_second_processing", exist_ok=True)


for root,subdir,files in os.walk("results_first_processing"):
    if len(files) > 100 :
        os.makedirs(os.path.join("results_second_processing",*root.split("/")[1:-1]), exist_ok=True)
        output_filename = os.path.join("results_second_processing",*root.split("/")[1:-1], root.split("/")[-1] + ".h5")

        with h5py.File(output_filename,"w") as outfile :

            n_runs = sum([file.endswith(".h5") for file in os.listdir(root)])
            crosses_00 = np.full((1000,n_runs),np.nan)
            theprod_01 = np.full((1000,n_runs),np.nan)
            r2s_00 = np.full((1000,n_runs),np.nan)
            r2s_vacancy = np.full((1000,n_runs),np.nan)
            theprod_00 = np.full((1000,n_runs),np.nan)
            theprod_11 = np.full((1000,n_runs),np.nan)
            r2s_11 = np.full((1000,n_runs),np.nan)
            crosses_11 = np.full((1000,n_runs),np.nan)

            i = 0
            for file in files:
                if file.endswith(".h5"):
                    input_filename = os.path.join(root,file)
                    with h5py.File(input_filename, "r") as f:
                        crosses_00[:, i] = f["cross_00"][:]
                        theprod_01[:, i] = f["theprod_01"][:]
                        r2s_00[:, i] = f["r2s_00"][:]
                        r2s_vacancy[:, i] = f["r2s_vacancy"][:]
                        theprod_00[:, i] = f["theprod_00"][:]
                        theprod_11[:, i] = f["theprod_11"][:]
                        r2s_11[:, i] = f["r2s_11"][:]
                        crosses_11[:, i] = f["cross_11"][:]

                        if i == 0 :
                            f.copy("input", outfile)

                        i += 1


            outfile.create_dataset("crosses_00", data=crosses_00.mean(axis=1))
            outfile.create_dataset("theprod_01", data=theprod_01.mean(axis=1))
            outfile.create_dataset("r2s_00", data=r2s_00.mean(axis=1))
            outfile.create_dataset("r2s_vacancy", data=r2s_vacancy.mean(axis=1))
            outfile.create_dataset("theprod_00", data=theprod_00.mean(axis=1))
            outfile.create_dataset("theprod_11", data=theprod_11.mean(axis=1))
            outfile.create_dataset("r2s_11", data=r2s_11.mean(axis=1))
            outfile.create_dataset("crosses_11", data=crosses_11.mean(axis=1))

            outfile.create_dataset("crosses_00_std", data=crosses_00.std(axis=1))
            outfile.create_dataset("theprod_01_std", data=theprod_01.std(axis=1))
            outfile.create_dataset("r2s_00_std", data=r2s_00.std(axis=1))
            outfile.create_dataset("r2s_vacancy_std", data=r2s_vacancy.std(axis=1))
            outfile.create_dataset("theprod_00_std", data=theprod_00.std(axis=1))
            outfile.create_dataset("theprod_11_std", data=theprod_11.std(axis=1))
            outfile.create_dataset("r2s_11_std", data=r2s_11.std(axis=1))
            outfile.create_dataset("crosses_11_std", data=crosses_11.std(axis=1))

            outfile.create_dataset("n_runs", data=n_runs)