#!/bin/bash

# List of elements to iterate over "Hapto" "general" 
elements=("101" "102" "105" "106" "107" "108"
"110" "111" "20" "27" "33" "36" "37" "38" "39"
"40" "41" "42" "43" "44" "45" "46" "47" "48" "49"
"50" "51" "52" "53" "54" "57" "58" "59"
"60" "61" "63" "66" "67" "68" "70" "71" 
"72" "73" "74" "75" "76" "80" "81" "82"
"83" "86" "87" "88" "89" "90" "94" "95" "98" "99")

# Number of times to repeat each element
num_repeats=50

# Loop over each element
for element in "${elements[@]}"; do
    echo "Processing element: $element"

    # Loop to submit multiple jobs for the current element
    for ((i=1; i<=$num_repeats; i++)); do
        # Define the Slurm job parameters
        job_name="${element}_job_$i"
        output_file="${element}_job_$i.out"
        error_file="${element}_job_$i.err"

        # Submit the Slurm job
        sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=$job_name
#SBATCH --output=$output_file
#SBATCH --error=$error_file
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=SiCore

echo "Running job for element: $element (Iteration: $i)"
# Add your command, job execution here
python3 run_estimation.py general "$i" "$element"
EOT

        echo "Submitted job: $job_name"
    done
done
