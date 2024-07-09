Analysis of Context Specific Metabolic Models
=============================================

This notebook can be downloaded from `this
link <https://github.com/KarthikDani/PHCCOProject/blob/main/gsmm/tutorials/analyse_csm.ipynb>`__

Define ``model_names`` with their associated paths in a ``dict``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    model_paths = {
        "epithelial_csm": "epithelial_csm.xml",
        "mesenchymal_csm": "mesenchymal_csm.xml"
    }

Import analysis module from ``gsmm``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``analyse_and_save_fluxes(..)`` is a function that carries all the
analysis and saves the data for further analysis to be utilised by
``visualisation`` module

.. code:: ipython3

    from gsmm.csm.analyse_csm import analyse_and_save_fluxes

.. code:: ipython3

    analyse_and_save_fluxes(model_paths)


.. parsed-literal::

    2024-07-10 01:15:54,415 - INFO - The current solver interface glpk doesn't support setting the optimality tolerance.


.. parsed-literal::

    Loading model: epithelial_csm from epithelial_csm.xml...


.. parsed-literal::

    2024-07-10 01:15:54,857 - INFO - Loaded model: epithelial_csm from epithelial_csm.xml
    2024-07-10 01:15:55,806 - INFO - The current solver interface glpk doesn't support setting the optimality tolerance.


.. parsed-literal::

    Successfully loaded model: epithelial_csm
    Loading model: mesenchymal_csm from mesenchymal_csm.xml...


.. parsed-literal::

    2024-07-10 01:15:56,280 - INFO - Loaded model: mesenchymal_csm from mesenchymal_csm.xml
    /Library/Frameworks/Python.framework/Versions/3.11/lib/python3.11/site-packages/gsmm/csm/analyse_csm.py:130: FutureWarning: A value is trying to be set on a copy of a DataFrame or Series through chained assignment using an inplace method.
    The behavior will change in pandas 3.0. This inplace method will never work because the intermediate object on which we are setting values always behaves as a copy.
    
    For example, when doing 'df[col].method(value, inplace=True)', try using 'df.method({col: value}, inplace=True)' or df[col] = df[col].method(value) instead, to perform the operation inplace on the original object.
    
    
      df['Flux'].fillna(0, inplace=True)
    2024-07-10 01:15:56,599 - WARNING - Metabolite damp_c not found in model epithelial_csm
    2024-07-10 01:15:56,608 - WARNING - Metabolite mag__hs_c not found in model epithelial_csm
    2024-07-10 01:15:56,609 - WARNING - Metabolite lpchol_hs_c not found in model epithelial_csm
    2024-07-10 01:15:56,609 - WARNING - Metabolite hdcea_c not found in model epithelial_csm
    2024-07-10 01:15:56,609 - WARNING - Metabolite hdca_c not found in model epithelial_csm
    2024-07-10 01:15:56,609 - WARNING - Metabolite ocdcea_c not found in model epithelial_csm
    2024-07-10 01:15:56,610 - WARNING - Metabolite ocdca_c not found in model epithelial_csm
    2024-07-10 01:15:56,610 - WARNING - Metabolite ptrc_c not found in model epithelial_csm
    2024-07-10 01:15:56,611 - WARNING - Metabolite spmd_c not found in model epithelial_csm
    2024-07-10 01:15:56,611 - WARNING - Metabolite sprm_c not found in model epithelial_csm
    2024-07-10 01:15:56,619 - WARNING - Metabolite Q10_c not found in model epithelial_csm
    2024-07-10 01:15:56,619 - WARNING - Metabolite paps_c not found in model epithelial_csm


.. parsed-literal::

    Successfully loaded model: mesenchymal_csm
    Extracting fluxes for model: epithelial_csm...
    Optimizing model: Recon3D...
    Optimization completed for model: Recon3D
    Extracting fluxes for model: mesenchymal_csm...
    Optimizing model: Recon3D...
    Optimization completed for model: Recon3D
    Flux extraction complete.
    Filtering fluxes with threshold: 1...
    Filtering complete. Number of reactions above threshold: 271
    Collecting sink fluxes for model: epithelial_csm...
    Warning: Metabolite damp_c not found in model epithelial_csm
    Warning: Metabolite mag__hs_c not found in model epithelial_csm
    Warning: Metabolite lpchol_hs_c not found in model epithelial_csm
    Warning: Metabolite hdcea_c not found in model epithelial_csm
    Warning: Metabolite hdca_c not found in model epithelial_csm
    Warning: Metabolite ocdcea_c not found in model epithelial_csm
    Warning: Metabolite ocdca_c not found in model epithelial_csm
    Warning: Metabolite ptrc_c not found in model epithelial_csm
    Warning: Metabolite spmd_c not found in model epithelial_csm
    Warning: Metabolite sprm_c not found in model epithelial_csm
    Warning: Metabolite Q10_c not found in model epithelial_csm
    Warning: Metabolite paps_c not found in model epithelial_csm
    Extracted sink fluxes for model: epithelial_csm
    Collecting sink fluxes for model: mesenchymal_csm...


.. parsed-literal::

    2024-07-10 01:15:56,760 - WARNING - Metabolite damp_c not found in model mesenchymal_csm
    2024-07-10 01:15:56,768 - WARNING - Metabolite mag__hs_c not found in model mesenchymal_csm
    2024-07-10 01:15:56,768 - WARNING - Metabolite lpchol_hs_c not found in model mesenchymal_csm
    2024-07-10 01:15:56,768 - WARNING - Metabolite hdcea_c not found in model mesenchymal_csm
    2024-07-10 01:15:56,776 - WARNING - Metabolite ocdcea_c not found in model mesenchymal_csm
    2024-07-10 01:15:56,777 - WARNING - Metabolite ocdca_c not found in model mesenchymal_csm
    2024-07-10 01:15:56,777 - WARNING - Metabolite spmd_c not found in model mesenchymal_csm
    2024-07-10 01:15:56,777 - WARNING - Metabolite sprm_c not found in model mesenchymal_csm
    2024-07-10 01:15:56,778 - WARNING - Metabolite Q10_c not found in model mesenchymal_csm
    2024-07-10 01:15:56,778 - WARNING - Metabolite paps_c not found in model mesenchymal_csm
    2024-07-10 01:15:56,782 - INFO - Data saved to flux_data.pkl and sink_flux_data.pkl


.. parsed-literal::

    Warning: Metabolite damp_c not found in model mesenchymal_csm
    Warning: Metabolite mag__hs_c not found in model mesenchymal_csm
    Warning: Metabolite lpchol_hs_c not found in model mesenchymal_csm
    Warning: Metabolite hdcea_c not found in model mesenchymal_csm
    Warning: Metabolite ocdcea_c not found in model mesenchymal_csm
    Warning: Metabolite ocdca_c not found in model mesenchymal_csm
    Warning: Metabolite spmd_c not found in model mesenchymal_csm
    Warning: Metabolite sprm_c not found in model mesenchymal_csm
    Warning: Metabolite Q10_c not found in model mesenchymal_csm
    Warning: Metabolite paps_c not found in model mesenchymal_csm
    Extracted sink fluxes for model: mesenchymal_csm
    Sink flux collection complete.
    Saving data to flux_data.pkl...
    Data saved successfully to flux_data.pkl and sink_flux_data.pkl


.. code:: ipython3

    from gsmm.csm.visualisation import plot_fluxes

Get Flux related plots across Context Specific Models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Clustermap for All reactions when taken into consideration for
   comparision across different context-specific models, as specified in
   ``model_paths`` above.
2. Correlation coefficients between each of the models when all the
   common reactions are considered.
3. Correlation for common sink reactions between each of the models.

.. code:: ipython3

    plot_fluxes('flux_data.pkl',
                'sink_flux_data.pkl',
                True)


.. parsed-literal::

    2024-07-10 01:16:04,878 - INFO - Loaded data from flux_data.pkl
    2024-07-10 01:16:04,880 - INFO - Loaded data from sink_flux_data.pkl
    2024-07-10 01:16:04,881 - INFO - Generating flux distribution clustermap...
    2024-07-10 01:16:05,316 - INFO - Flux distribution clustermap saved as flux_distribution_clustermap.png



.. image:: analyse_csm_files/analyse_csm_9_1.png


.. parsed-literal::

    2024-07-10 01:16:05,474 - INFO - Generating flux correlation heatmap...
    2024-07-10 01:16:05,733 - INFO - Flux correlation heatmap saved as flux_correlation_heatmap.png



.. image:: analyse_csm_files/analyse_csm_9_3.png


.. parsed-literal::

    2024-07-10 01:16:05,810 - INFO - Generating sink fluxes heatmap...
    2024-07-10 01:16:06,039 - INFO - Sink fluxes heatmap saved as sink_fluxes_heatmap.png
    2024-07-10 01:16:06,039 - INFO - Generating sink flux correlation heatmap...
    2024-07-10 01:16:06,161 - INFO - Sink flux correlation heatmap saved as sink_flux_correlation_heatmap.png



.. image:: analyse_csm_files/analyse_csm_9_5.png


Similarly, two or more models can be compared to get the relevant plot
for significant observations!
