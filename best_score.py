import os
import json

def find_file_with_lowest_best_score(directory):
    lowest_score = float('inf')  # Set initial value to infinity
    file_with_lowest_score = None

    # Iterate over all files in the directory
    for filename in os.listdir(directory):
        #print(filename)
        if filename.endswith(".json"):  # Check if the file is a JSON file
            file_path = os.path.join(directory, filename)
            #print(filename)
            if not filename.startswith("2024_10_"):continue
            try:
                # Open and load the JSON file as a dictionary
                with open(file_path, 'r') as f:
                    data = json.load(f)
                if "Hkt_init" in data['best_parameters']: continue
                # Check if the dictionary contains the "best_score" key
                if 'best_score' in data:
                    best_score = data['best_score']

                    # Compare and find the file with the lowest best_score
                    if best_score < lowest_score:
                        lowest_score = best_score
                        file_with_lowest_score = filename
            except Exception as e:
                print(f"Error reading {file_path}: {e}")

    # Return the file with the lowest best_score and the score itself
    return file_with_lowest_score, lowest_score
def main():
    # Example usage:
    for model in ['general','Hapto','immune']:
        directory_path = f"Estimation/{model}/pth"
        print(directory_path)
        file, score = find_file_with_lowest_best_score(directory_path)
        if file:
            print(f"The file with the lowest best_score is {file} with a score of {score} in {model}.")
        else:
            print("No file with a 'best_score' found.")

if __name__=='__main__':
    main()
