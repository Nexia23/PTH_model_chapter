#!/bin/bash

# List of elements to iterate over
elements=("Hapto" "immune" "general")

# Number of times to repeat each element
num_repeats=10

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
# Add your command or job execution here
# Example: Assuming you want to run a Python script with the element as an argument
python3 run_estimation.py "$element" "$i"
EOT

        echo "Submitted job: $job_name"
    done
done
