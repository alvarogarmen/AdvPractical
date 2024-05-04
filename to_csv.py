import re
import csv
import sys

# Check if correct number of command-line arguments are provided
if len(sys.argv) != 3:
    print("Usage: python script_name.py input_file output_file")
    sys.exit(1)

# Get input and output file names from command-line arguments
input_file_name = sys.argv[1]
output_file_name = sys.argv[2]

# Define regular expressions for extracting information
hypergraph_name_pattern = re.compile(r'name: "(.*?)"')
node_count_pattern = re.compile(r'node_count: (\d+)')
nanos_pattern = re.compile(r'nanos: (\d+)')
weight_pattern = re.compile(r'weight: (\d+)')

# Initialize flag to indicate if currently processing a hypergraph block
processing_hypergraph = False

# Open input file in read mode
with open(input_file_name, 'r') as input_file:
    # Open output file in write mode
    with open(output_file_name, 'w', newline='') as output_file:
        # Create CSV writer
        writer = csv.writer(output_file)
        
        # Write header to CSV file
        writer.writerow(['Hypergraph Name', 'Node Count', 'Algorithm Duration (nanos)', 'Weight'])
        
        # Initialize variables to store information
        hypergraph_name = None
        node_count = None
        nanos = None
        weight = None
        
        # Read input file line by line
        for line in input_file:
            # Check if line starts with 'results {'
            if line.strip() == 'results {':
                # Reset variable for the new block
                hypergraph_name = None
                processing_hypergraph = False
            elif line.strip().startswith('hypergraph {'):
                processing_hypergraph = True
            elif processing_hypergraph:
                # Use regular expression to extract hypergraph name
                match = hypergraph_name_pattern.search(line)
                if match:
                    hypergraph_name = match.group(1)
                    print("Hypergraph Name:", hypergraph_name)
                    # Reset flag to stop further processing within hypergraph block
                    processing_hypergraph = False
            else:
                # Use regular expressions to extract information within the results block
                match = node_count_pattern.search(line)
                if match:
                    node_count = match.group(1)
                    print("Node Count:", node_count)
                
                match = nanos_pattern.search(line)
                if match:
                    nanos = match.group(1)
                    print("Algorithm Duration (nanos):", nanos)
                
                match = weight_pattern.search(line)
                if match:
                    weight = match.group(1)
                    print("Weight:", weight)
        
                # If all fields are found, write to CSV
                if hypergraph_name and node_count and nanos and weight:
                    writer.writerow([hypergraph_name, node_count, nanos, weight])
                    hypergraph_name = None
                    node_count = None
                    nanos = None
                    weight = None