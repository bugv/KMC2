import numpy as np
import h5py 
import os 
import copy

os.makedirs("results_second_processing", exist_ok=True)


for root,subdir,files in os.walk("results_first_processing"):
    if len(files) > 100 :
        print("Processing", root)
        print("hi")
        os.makedirs(os.path.join("results_second_processing",*root.split("/")[1:-1]), exist_ok=True)
        output_filename = os.path.join("results_second_processing",*root.split("/")[1:-1], root.split("/")[-1] + ".h5")

        with h5py.File(output_filename,"w") as outfile :

            n_runs = sum([file.endswith(".h5") for file in os.listdir(root)])
            crosses_AA = np.full((1000,n_runs),np.nan)
            theprod_AA = np.full((1000,n_runs),np.nan)
            r2s_AA = np.full((1000,n_runs),np.nan)
            r2s_vacancy = np.full((1000,n_runs),np.nan)
            theprod_BB = np.full((1000,n_runs),np.nan)
            r2s_BB = np.full((1000,n_runs),np.nan)
            crosses_BB = np.full((1000,n_runs),np.nan)
            theprod_AB = np.full((1000,n_runs),np.nan)
            time_collector = np.full((1000,n_runs),np.nan)


            i = 0
            for file in files:
                if file.endswith(".h5"):
                    input_filename = os.path.join(root,file)
                    with h5py.File(input_filename, "r") as f:
                        if f["input"]["atom_key"]["Al"][()] == 0 and f["input"]["atom_key"]["Fe"][()] == 1 :
                            crosses_AA[:, i] = f["cross_00"][:]
                            theprod_AA[:, i] = f["theprod_00"][:]
                            r2s_AA[:, i] = f["r2s_00"][:]

                            r2s_vacancy[:, i] = f["r2s_vacancy"][:]
                            theprod_AB[:, i] = f["theprod_01"][:]

                            theprod_BB[:, i] = f["theprod_11"][:]
                            r2s_BB[:, i] = f["r2s_11"][:]
                            crosses_BB[:, i] = f["cross_11"][:]

                        elif f["input"]["atom_key"]["Fe"][()] == 0 and f["input"]["atom_key"]["Al"][()] == 1 :

                            crosses_BB[:, i] = f["cross_00"][:]
                            theprod_BB[:, i] = f["theprod_00"][:]
                            r2s_BB[:, i] = f["r2s_00"][:]

                            r2s_vacancy[:, i] = f["r2s_vacancy"][:]
                            theprod_AB[:, i] = f["theprod_01"][:]

                            theprod_AA[:, i] = f["theprod_11"][:]
                            r2s_AA[:, i] = f["r2s_11"][:]
                            crosses_AA[:, i] = f["cross_11"][:]
                        else : 
                            print(f"Error: atom_key not recognized in {input_filename}")
                        
                        time_collector[:, i] = f["time_collector"][:]

                        if i == 0 :
                            f.copy("input", outfile)

                        i += 1


            outfile.create_dataset("crosses_AA", data=crosses_AA.mean(axis=1))
            outfile.create_dataset("theprod_AB", data=theprod_AB.mean(axis=1))
            outfile.create_dataset("r2s_AA", data=r2s_AA.mean(axis=1))
            outfile.create_dataset("r2s_vacancy", data=r2s_vacancy.mean(axis=1))
            outfile.create_dataset("theprod_AA", data=theprod_AA.mean(axis=1))
            outfile.create_dataset("theprod_BB", data=theprod_BB.mean(axis=1))
            outfile.create_dataset("r2s_BB", data=r2s_BB.mean(axis=1))
            outfile.create_dataset("crosses_BB", data=crosses_BB.mean(axis=1))
            outfile.create_dataset("time_collector", data=time_collector.mean(axis=1))

            outfile.create_dataset("crosses_AA_std", data=crosses_AA.std(axis=1))
            outfile.create_dataset("theprod_AB_std", data=theprod_AB.std(axis=1))
            outfile.create_dataset("r2s_AA_std", data=r2s_AA.std(axis=1))
            outfile.create_dataset("r2s_vacancy_std", data=r2s_vacancy.std(axis=1))
            outfile.create_dataset("theprod_AA_std", data=theprod_AA.std(axis=1))
            outfile.create_dataset("theprod_BB_std", data=theprod_BB.std(axis=1))
            outfile.create_dataset("r2s_BB_std", data=r2s_BB.std(axis=1))
            outfile.create_dataset("crosses_BB_std", data=crosses_BB.std(axis=1))
            outfile.create_dataset("time_collector_std", data=time_collector.std(axis=1))

            outfile.create_dataset("n_runs", data=n_runs)