Building Context Specific Metabolic Models
==========================================

This notebook can be downloaded from `this
link <https://github.com/KarthikDani/PHCCOProject/blob/main/gsmm/tutorials/build_csm.ipynb>`__

-  ``run_model_reconstruction`` from ``gsmm.csm.build_csm``: Used for
   running the model reconstruction pipeline.
-  ``cobra``: A package for constraint-based reconstruction and
   analysis.

.. code:: ipython3

    from gsmm.csm.build_csm import run_model_reconstruction
    
    import cobra
    import pandas as pd

Load expression data
~~~~~~~~~~~~~~~~~~~~

For example, let’s an EMT expression data to follow along this tutorial
- ``emt_expression_data_path``: Path to the CSV file containing the EMT
expression data. - ``emt_expression_data``: DataFrame holding the loaded
EMT expression data.

.. code:: ipython3

    emt_expression_data_path = "../../Data/InputData/EMT_FINAL_DATA.csv"
    emt_expression_data = pd.read_csv(emt_expression_data_path)
    emt_expression_data.head()




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>Gene_ID</th>
          <th>Mesenchymal</th>
          <th>Epithelial</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>2978_AT1</td>
          <td>1.590799</td>
          <td>-7.624877</td>
        </tr>
        <tr>
          <th>1</th>
          <td>1571_AT1</td>
          <td>6.046029</td>
          <td>3.480441</td>
        </tr>
        <tr>
          <th>2</th>
          <td>1549_AT1</td>
          <td>2.358338</td>
          <td>3.101581</td>
        </tr>
        <tr>
          <th>3</th>
          <td>1548_AT1</td>
          <td>2.358338</td>
          <td>3.101581</td>
        </tr>
        <tr>
          <th>4</th>
          <td>949_AT1</td>
          <td>8.881953</td>
          <td>8.613039</td>
        </tr>
      </tbody>
    </table>
    </div>



Load Recon3D model from the web
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  ``recon_model``: The Recon3D model loaded from a web source using
   COBRApy.

This step loads the Recon3D model, a comprehensive genome-scale
metabolic model. Recon3D integrates metabolic data from various human
tissues and is widely used for studying human metabolism.

Structure of a COBRApy SBML Model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. figure:: cobrapy_model_structure.png
   :alt: Cobrapy Model Structure

   Cobrapy Model Structure

In COBRApy, a SBML model is structured as follows:

-  ``Model``: The core object representing the metabolic model. It
   contains several attributes and methods to manipulate and analyze the
   model.

   -  **Attributes**:

      -  ``reactions``: A list of ``Reaction`` objects representing the
         biochemical reactions in the model.
      -  ``metabolites``: A list of ``Metabolite`` objects representing
         the metabolites in the model.
      -  ``genes``: A list of ``Gene`` objects representing the genes
         associated with the reactions.
      -  ``objective``: The objective function of the model, typically
         used in flux balance analysis (FBA).
      -  ``id``: A unique identifier for the model.
      -  ``name``: A descriptive name for the model.

   -  **Methods**:

      -  ``optimize()``: Performs flux balance analysis to optimize the
         objective function.
      -  ``summary()``: Provides a summary of the model, including the
         number of reactions, metabolites, and genes.

-  ``Reaction``: Represents a biochemical reaction in the model.

   -  **Attributes**:

      -  ``id``: A unique identifier for the reaction.
      -  ``name``: A descriptive name for the reaction.
      -  ``metabolites``: A dictionary of ``Metabolite`` objects and
         their stoichiometric coefficients in the reaction.
      -  ``lower_bound``: The lower bound of the reaction flux.
      -  ``upper_bound``: The upper bound of the reaction flux.

   -  **Methods**:

      -  ``add_metabolites(metabolites)``: Adds metabolites to the
         reaction.
      -  ``remove_metabolites(metabolites)``: Removes metabolites from
         the reaction.

-  ``Metabolite``: Represents a metabolite in the model.

   -  **Attributes**:

      -  ``id``: A unique identifier for the metabolite.
      -  ``name``: A descriptive name for the metabolite.
      -  ``formula``: The chemical formula of the metabolite.
      -  ``compartment``: The compartment where the metabolite is
         located.

-  ``Gene``: Represents a gene in the model.

   -  **Attributes**:

      -  ``id``: A unique identifier for the gene.
      -  ``name``: A descriptive name for the gene.
      -  ``reactions``: A list of ``Reaction`` objects associated with
         the gene.

Now we load ``Recon3D`` model as the parent model from which the base model will be created following ``gsmm`` library’s functions, to further reconstruct context specific models be done.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    recon_model = cobra.io.web.load_model(model_id="Recon3D")
    recon_model


.. parsed-literal::

    2024-07-09 22:15:26,787 - INFO - The current solver interface glpk doesn't support setting the optimality tolerance.




.. raw:: html

    
    <table>
        <tr>
            <td><strong>Name</strong></td>
            <td>Recon3D</td>
        </tr><tr>
            <td><strong>Memory address</strong></td>
            <td>106b88ed0</td>
        </tr><tr>
            <td><strong>Number of metabolites</strong></td>
            <td>5835</td>
        </tr><tr>
            <td><strong>Number of reactions</strong></td>
            <td>10600</td>
        </tr><tr>
            <td><strong>Number of genes</strong></td>
            <td>2248</td>
        </tr><tr>
            <td><strong>Number of groups</strong></td>
            <td>0</td>
        </tr><tr>
            <td><strong>Objective expression</strong></td>
            <td>1.0*BIOMASS_maintenance - 1.0*BIOMASS_maintenance_reverse_5b3f9</td>
        </tr><tr>
            <td><strong>Compartments</strong></td>
            <td>cytosol, lysosome, mitochondria, endoplasmic reticulum, extracellular space, peroxisome/glyoxysome, nucleus, golgi apparatus, inner mitochondrial compartment</td>
        </tr>
      </table>



We save ``Recon3D`` as ``.xml`` file.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    cobra.io.write_sbml_model(recon_model, "recon_model.xml")

Reconstruction algorithm in ``gsmm`` needs, 1. ``parent_model_path``:
``Recon3D`` in this case. 2. ``base_model_path``: All the unneccessary
genes, reactions and metabolites are removed, then the model is saved as
``base_model``. 3. ``gene_id_column``: Expression data in our case has
``gene ids`` in the column called ``Gene_ID``. We need to tell the
algorithm that this is our ``gene id`` column.

.. code:: ipython3

    parent_model_path = "recon_model.xml"
    base_model_path = "recon_model.xml"
    gene_id_column = "Gene_ID"

Running Model Reconstruction to get Context-Specific Model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``run_model_reconstruction(..)`` takes care of the interface to run the
main pipeline for reconstructing an optimized metabolic model from the
provided expression data. It handles exceptions and prints error
messages if reconstruction fails, returning None in case of errors.

.. code:: ipython3

    epithelial_csm = run_model_reconstruction(model_path=parent_model_path,
                                              base_model_path=base_model_path,
                                              data_path=emt_expression_data_path,
                                              gene_id_column=gene_id_column,
                                              scores_column="Epithelial",
                                              )
    epithelial_csm


.. parsed-literal::

    2024-07-09 23:07:13,201 - INFO - The current solver interface glpk doesn't support setting the optimality tolerance.


.. parsed-literal::

    Cobra Configuration set to GLPK solver!
    Loading SBML model from recon_model.xml...


.. parsed-literal::

    2024-07-09 23:07:16,085 - INFO - The current solver interface glpk doesn't support setting the optimality tolerance.


.. parsed-literal::

    SBML model loaded.
    Loading expression data from ../../Data/InputData/EMT_FINAL_DATA.csv...
    Expression data loaded.
    Extracting genes from Gene_ID column...
    Number of Genes extracted: 1786
    Filtering model by genes...


.. parsed-literal::

    2024-07-09 23:07:17,796 - INFO - The current solver interface glpk doesn't support setting the optimality tolerance.
    2024-07-09 23:07:18,691 - INFO - The current solver interface glpk doesn't support setting the optimality tolerance.


.. parsed-literal::

    Model filtered.
    Normalizing expression data...
    Expression data has been normalized.
    Assigning confidence levels to reactions giving highest confidence levels to Biomass reactions...


.. parsed-literal::

    2024-07-09 23:07:25,589 - INFO - The current solver interface glpk doesn't support setting the optimality tolerance.


.. parsed-literal::

    Reaction Confidences:                 Reaction_ID  Confidence_Level
    0           24_25DHVITD3tm                -1
    1                25HVITD3t                -1
    2                    COAtl                -1
    3       EX_5adtststerone_e                -1
    4      EX_5adtststerones_e                -1
    ...                    ...               ...
    10436         ACMPGLUTTRsc                -1
    10437             FVSCOAhc                -1
    10438             MDZGLChr                -1
    10439             TMACMPhr                -1
    10440           CYSACMPitr                -1
    
    [10441 rows x 2 columns]
    Reaction confidence levels saved as reaction_Epithelial_confidence_levels.csv.
    Initializing and building CORDA model...


.. parsed-literal::

    2024-07-09 23:37:08,413 - INFO - The current solver interface glpk doesn't support setting the optimality tolerance.


.. parsed-literal::

    CORDA model optimization completed.
    Optimized CORDA model:
    build status: reconstruction complete
    Inc. reactions: 1539/10441
     - unclear: 213/1071
     - exclude: 738/6995
     - low and medium: 360/2145
     - high: 228/230
    
    Reconstruction completed from the parent model Recon3D




.. raw:: html

    
    <table>
        <tr>
            <td><strong>Name</strong></td>
            <td>Recon3D</td>
        </tr><tr>
            <td><strong>Memory address</strong></td>
            <td>3317be8d0</td>
        </tr><tr>
            <td><strong>Number of metabolites</strong></td>
            <td>1198</td>
        </tr><tr>
            <td><strong>Number of reactions</strong></td>
            <td>1539</td>
        </tr><tr>
            <td><strong>Number of genes</strong></td>
            <td>683</td>
        </tr><tr>
            <td><strong>Number of groups</strong></td>
            <td>0</td>
        </tr><tr>
            <td><strong>Objective expression</strong></td>
            <td>1.0*BIOMASS_maintenance - 1.0*BIOMASS_maintenance_reverse_5b3f9</td>
        </tr><tr>
            <td><strong>Compartments</strong></td>
            <td>cytosol, mitochondria, peroxisome/glyoxysome, extracellular space, golgi apparatus, lysosome, endoplasmic reticulum, nucleus</td>
        </tr>
      </table>



Similarly we do so to get ``mesenchymal`` context specific metabolic
model.

.. code:: ipython3

    mesenchymal_csm = run_model_reconstruction(model_path=parent_model_path,
                                              base_model_path=base_model_path,
                                              data_path=emt_expression_data_path,
                                              gene_id_column=gene_id_column,
                                              scores_column="Mesenchymal",
                                              )
    mesenchymal_csm


.. parsed-literal::

    2024-07-09 23:37:41,993 - INFO - The current solver interface glpk doesn't support setting the optimality tolerance.


.. parsed-literal::

    Cobra Configuration set to GLPK solver!
    Loading SBML model from recon_model.xml...


.. parsed-literal::

    2024-07-09 23:37:45,448 - INFO - The current solver interface glpk doesn't support setting the optimality tolerance.


.. parsed-literal::

    SBML model loaded.
    Loading expression data from ../../Data/InputData/EMT_FINAL_DATA.csv...
    Expression data loaded.
    Extracting genes from Gene_ID column...
    Number of Genes extracted: 1786
    Filtering model by genes...


.. parsed-literal::

    2024-07-09 23:37:46,399 - INFO - The current solver interface glpk doesn't support setting the optimality tolerance.
    2024-07-09 23:37:47,240 - INFO - The current solver interface glpk doesn't support setting the optimality tolerance.


.. parsed-literal::

    Model filtered.
    Normalizing expression data...
    Expression data has been normalized.
    Assigning confidence levels to reactions giving highest confidence levels to Biomass reactions...


.. parsed-literal::

    2024-07-09 23:37:54,345 - INFO - The current solver interface glpk doesn't support setting the optimality tolerance.


.. parsed-literal::

    Reaction Confidences:                 Reaction_ID  Confidence_Level
    0           24_25DHVITD3tm                -1
    1                25HVITD3t                -1
    2                    COAtl                -1
    3       EX_5adtststerone_e                -1
    4      EX_5adtststerones_e                -1
    ...                    ...               ...
    10436         ACMPGLUTTRsc                -1
    10437             FVSCOAhc                -1
    10438             MDZGLChr                -1
    10439             TMACMPhr                -1
    10440           CYSACMPitr                -1
    
    [10441 rows x 2 columns]
    Reaction confidence levels saved as reaction_Mesenchymal_confidence_levels.csv.
    Initializing and building CORDA model...


.. parsed-literal::

    2024-07-10 00:06:41,053 - INFO - The current solver interface glpk doesn't support setting the optimality tolerance.


.. parsed-literal::

    CORDA model optimization completed.
    Optimized CORDA model:
    build status: reconstruction complete
    Inc. reactions: 1518/10441
     - unclear: 221/1523
     - exclude: 713/6606
     - low and medium: 385/2107
     - high: 199/205
    
    Reconstruction completed from the parent model Recon3D




.. raw:: html

    
    <table>
        <tr>
            <td><strong>Name</strong></td>
            <td>Recon3D</td>
        </tr><tr>
            <td><strong>Memory address</strong></td>
            <td>32ddb6490</td>
        </tr><tr>
            <td><strong>Number of metabolites</strong></td>
            <td>1200</td>
        </tr><tr>
            <td><strong>Number of reactions</strong></td>
            <td>1518</td>
        </tr><tr>
            <td><strong>Number of genes</strong></td>
            <td>684</td>
        </tr><tr>
            <td><strong>Number of groups</strong></td>
            <td>0</td>
        </tr><tr>
            <td><strong>Objective expression</strong></td>
            <td>1.0*BIOMASS_maintenance - 1.0*BIOMASS_maintenance_reverse_5b3f9</td>
        </tr><tr>
            <td><strong>Compartments</strong></td>
            <td>cytosol, mitochondria, peroxisome/glyoxysome, extracellular space, lysosome, endoplasmic reticulum, nucleus</td>
        </tr>
      </table>



Save the CSMs for further analysis and plots
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    cobra.io.write_sbml_model(epithelial_csm, "epithelial_csm.xml")
    cobra.io.write_sbml_model(mesenchymal_csm, "mesenchymal_csm.xml")
