import os
import subprocess

# Define the path to your .ipynb file
notebook_path = "csm_analysis_.ipynb"
base_filename = os.path.basename(notebook_path).split(".")[0]
markdown_path = os.path.join(os.path.dirname(notebook_path), base_filename + ".md")
pdf_output_path = os.path.join(os.path.dirname(notebook_path), base_filename + ".pdf")

metadata = {
    "title": "Analysing Flux Chnges, their Correlation and Nutrient Depletion with Amino acids as priority.",
    "author": "Karthik",
    "date": "2024-07-07",
    "abstract": "Context-specific metabolic models (CSMs) are crucial for representing metabolic processes in particular biological conditions. This report outlines the methodology for constructing CSMs tailored to epithelial and mesenchymal states using the Open source CORDA algorithms, integrated with RNA-Seq data. We describe the steps taken to refine generic metabolic models, ensuring that they reflect the metabolic activities specific to the studied contexts. This approach provides a robust framework for building detailed and condition-specific models for further metabolic analysis.",
    #"subtitle": "Process of constructing CSMs for epithelial and mesenchymal states using the CORDA algorithm.",
    "keywords": "Genome Scale Metabolic Modeling",
    "thanks": "*",
    "email": "karthikdani14@gmail.com",
    "institute": "BMS College of Engineering, Bangalore"
}

# try:
#     subprocess.run(['jupyter', 'nbconvert', '--to', 'markdown', #'--TagRemovePreprocessor.enabled=True', '--TagRemovePreprocessor.remove_cell_tags="exclude-output', 
#                     notebook_path], check=True)
#     print(f"Converted notebook to Markdown: {markdown_path}")
    
# except subprocess.CalledProcessError as e:
#     print(f"Error during Markdown conversion: {e}")
#     exit(1)

pandoc_command = ["pandoc", notebook_path, "-o", pdf_output_path]

for key, value in metadata.items():
    pandoc_command.append(f"--metadata={key}={value}")

# Convert the Markdown file to PDF using pandoc
try:
    subprocess.run(pandoc_command, check=True)
    print(f"Successfully converted Markdown to PDF: {pdf_output_path}")
except subprocess.CalledProcessError as e:
    print(f"Error during PDF conversion: {e}")
    print(f"Error during PDF conversion: {e}")

# Clean up intermediate files if needed
# os.remove(markdown_path)
