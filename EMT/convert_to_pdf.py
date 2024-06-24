import nbconvert as nb
import os
import subprocess

# Define the path to your .ipynb file
#notebook_path = "/Users/karthik/Desktop/PHCCO IISc Internship/EMT/integrate_expression.ipynb"
notebook_path = "/Users/karthik/Desktop/PHCCO IISc Internship/EMT/fasting_simulation.ipynb"

input_file = notebook_path
output_file = notebook_path.split("/")[-1].split(".")[0] + ".pdf"

#title = 'Integrating Fasting expression data to the EMT SBML model'
title = 'Comparative Metabolic Analysis of Fasting Effects on an EMT Model'
authors = ["Karthik"]
date = "2024-06-23"
#abstract = "Contains integrated fasting-induced gene expression data with the epithelial-mesenchymal transition (EMT) signaling networks using computational methods. The EMT model, derived from SBML pathway maps and enhanced with COBRApy, was refined to incorporate gene regulatory relationships and metabolic constraints. Focusing on the impact of fasting on EMT metastasis, we identified key genes and reactions influenced by differential gene expression profiles. Through relaxation techniques and flux variability analysis, we simulated the metabolic responses under fasting conditions, revealing potential regulatory nodes and metabolic pathways affected during EMT. Our approach provides insights into how fasting-induced gene expression changes may modulate EMT dynamics, offering a computational framework for studying metabolic adaptations in cancer progression."
abstract = "This study employs metabolic modeling to investigate the impact of fasting on an Epithelial-Mesenchymal Transition model. Using constraint-based modeling techniques, we compare the metabolic flux distributions between an original model and a fasting integrated version for fasting conditions. Flux Variability Analysis is conducted to assess the range of feasible flux values for each reaction in both models, highlighting reactions significantly altered by fasting. This computational approach provides insights into metabolic adaptations during EMT under fasting conditions, potentially uncovering metabolic targets."

keywords = "Genome Scale Metabolic Modeling"

#subtitle = "Integrating Fasting-Induced Gene Expressions with EMT Signaling Networks"
subtitle = "Exploring Metabolic Adaptations in EMT Models under Fasting Conditions."

thanks = "Thanks"
affiliation = ["Affiliation 1"]
email = ["karthikdani14@gmail.com"]
institute = ["BMS College of Engineering, Bangalore"]

metadata = [
    f'--metadata=title={title}',
    f'--metadata=author={", ".join(authors)}',
    f'--metadata=date={date}',
    f'--metadata=abstract={abstract}',
    #f'--metadata=keywords={keywords}',
    f'--metadata=subtitle={subtitle}',
    #f'--metadata=thanks={thanks}',
    #f'--metadata=affiliation={", ".join(affiliation)}',
    #f'--metadata=email={", ".join(email)}',
    #f'--metadata=institute={", ".join(institute)}'
]

pandoc_command = ["pandoc", input_file, "-o", output_file] + metadata

# for subcommand in pandoc_command:
#     print(subcommand, end=" ")

try:
    subprocess.run('jupyter nbconvert --to markdown ' + notebook_path)
    subprocess.run(pandoc_command, check=True)
    print(f"Successfully converted {input_file} to {output_file}")
except subprocess.CalledProcessError as e:
    print(f"Error during conversion: {e}")

# os.system(pandoc_command)