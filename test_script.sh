input_filename="inputs/composition_Al_0.50_Fe_0.50.json"
tmp_raw_kmc_output_file="tmp/results.json"
output_filename="results_first_processing/composition_Al_0.50_Fe_0.50/test.h5"

python main.py -binary $input_filename
python analysis_v2.py $tmp_raw_kmc_output_file $output_filename

