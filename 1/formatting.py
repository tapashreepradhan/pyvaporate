def format_and_write_file(input_file, output_file):
    """
    Reads a text file, reformats its numerical content while preserving the first and last lines,
    and writes the formatted data to a new file.

    Parameters:
        input_file (str): Path to the input file.
        output_file (str): Path to the output file.
    """
    with open(input_file, 'r') as f:
        lines = f.readlines()

    if len(lines) < 3:
        print("Error: File must have at least three lines.")
        return

    formatted_lines = []  # Initialize list for formatted content

    for line in lines[1:-1]:  # Exclude the first and last lines
        parts = line.strip().split()
        if len(parts) == 5:  # Ensure the line has exactly 5 elements
            formatted_line = f"{float(parts[0]):.10e}\t{float(parts[1]):.10e}\t{float(parts[2]):.10e}\t{parts[3]}\t{parts[4]}"
            formatted_lines.append(formatted_line)

    with open(output_file, 'w') as f:
        f.write(lines[0].strip() + "\n")  # Write the first line
        f.write("\n".join(formatted_lines) + "\n")  # Write formatted lines
        f.write(lines[-1].strip() + "\n")  # Write the last line

# Usage:
input_file = "./mesh.txt"   # Change this to your actual input file path
output_file = "./f_mesh.txt"  # Change this to your desired output file path

format_and_write_file(input_file, output_file)

print(f"Formatted data has been written to: {output_file}")