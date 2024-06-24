import os
import subprocess

# Define the path to your .ipynb file
notebook_path = "/Users/karthik/Desktop/PHCCO IISc Internship/EMT/integrate_expression.ipynb"
base_filename = os.path.basename(notebook_path).split(".")[0]
markdown_path = os.path.join(os.path.dirname(notebook_path), base_filename + ".md")
pdf_output_path = os.path.join(os.path.dirname(notebook_path), base_filename + ".pdf")

# Metadata for the PDF
metadata_comp = {
    "title": "Comparative Metabolic Analysis of Fasting Effects on an EMT Model",
    "author": "Karthik",
    "date": "2024-06-23",
    "abstract": ("This study employs metabolic modeling to investigate the impact of fasting on an "
                 "Epithelial-Mesenchymal Transition model. Using constraint-based modeling techniques, "
                 "we compare the metabolic flux distributions between an original model and a fasting "
                 "integrated version for fasting conditions. Flux Variability Analysis is conducted to "
                 "assess the range of feasible flux values for each reaction in both models, highlighting "
                 "reactions significantly altered by fasting. This computational approach provides insights "
                 "into metabolic adaptations during EMT under fasting conditions, potentially uncovering "
                 "metabolic targets."),
    "subtitle": "Exploring Metabolic Adaptations in EMT Models under Fasting Conditions.",
    "keywords": "Genome Scale Metabolic Modeling",
    "thanks": "Thanks",
    "email": "karthikdani14@gmail.com",
    "institute": "BMS College of Engineering, Bangalore"
}

metadata = {
    "title": "Integrating Fasting expression data to the EMT SBML model",
    "author": "Karthik",
    "date": "2024-06-23",
    "abstract": ("Contains integrated fasting-induced gene expression data with the epithelial-mesenchymal transition (EMT) signaling networks using computational methods. The EMT model, derived from SBML pathway maps and enhanced with COBRApy, was refined to incorporate gene regulatory relationships and metabolic constraints. Focusing on the impact of fasting on EMT metastasis, we identified key genes and reactions influenced by differential gene expression profiles. Through relaxation techniques and flux variability analysis, we simulated the metabolic responses under fasting conditions, revealing potential regulatory nodes and metabolic pathways affected during EMT. Our approach provides insights into how fasting-induced gene expression changes may modulate EMT dynamics, offering a computational framework for studying metabolic adaptations in cancer progression."),
    "subtitle": "Integrating Fasting-Induced Gene Expressions with EMT Signaling Networks",
    "keywords": "Genome Scale Metabolic Modeling",
    "thanks": "*",
    "email": "karthikdani14@gmail.com",
    "institute": "BMS College of Engineering, Bangalore"
}

# Convert the notebook to Markdown
try:
    subprocess.run(['jupyter', 'nbconvert', '--to', 'markdown', '--TagRemovePreprocessor.enabled=True', '--TagRemovePreprocessor.remove_cell_tags="exclude-output', notebook_path], check=True)
    print(f"Converted notebook to Markdown: {markdown_path}")
    
except subprocess.CalledProcessError as e:
    print(f"Error during Markdown conversion: {e}")
    exit(1)

# Build the pandoc command with metadata
pandoc_command = ["pandoc", markdown_path, "-o", pdf_output_path]

for key, value in metadata.items():
    pandoc_command.append(f"--metadata={key}={value}")

# Convert the Markdown file to PDF using pandoc
try:
    subprocess.run(pandoc_command, check=True)
    print(f"Successfully converted Markdown to PDF: {pdf_output_path}")
except subprocess.CalledProcessError as e:
    print(f"Error during PDF conversion: {e}")

# Clean up intermediate files if needed
# os.remove(markdown_path)
