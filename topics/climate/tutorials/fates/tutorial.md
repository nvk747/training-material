---
layout: tutorial_hands_on
title: Functionally Assembled Terrestrial Ecosystem Simulator (FATES)
enable: false
zenodo_link: 'https://doi.org/10.5281/zenodo.4108341'
requirements:
  -
    type: "internal"
    topic_name: climate
    tutorials:
        - panoply

questions:
- How to run CLM-FATES with the CLM-FATES Galaxy tool?
- How to upload input data for running CLM-FATES?
- How to customize your runs?
- How to analyze your model outputs?
- How to create a workflow?
- How to share your workflow?
objectives:
- Setting up a CLM-FATES case.
- Customizing your run.
- Interactive visualization with Panoply.
- Automating your analyzes and visualisations of your CLM-FATES case.
- Creating multi-case scenarios.
- Composing, executing and publishing CML-FATES workflow.
time_estimation: 4H
key_points:
- CLM-FATES
- Quick visualization of your results with Panoply
- Create multi-case simulations with a Galaxy workflow
contributors:
- annefou

---


# Introduction
{:.no_toc}


The practical aims at familiarizing you with running CLM-FATES in Galaxy and analyzing the model results.
It will also teach you on how to create Galaxy workflow for your CLM-FATES simulations to make your research fully reproducible.

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

> ### {% icon comment %} Background
>
> FATES is the “Functionally Assembled Terrestrial Ecosystem Simulator”.
> FATES needs what we call a "Host Land Model" (HLM) to run and in this tutorial
> we will be using the [Community Land Model](http://www.cesm.ucar.edu/models/clm/)
> of the [Community Terrestrial Systems Model](https://github.com/ESCOMP/CTSM) (CLM-CTSM).
> FATES was derived from the CLM Ecosystem Demography model (CLM(ED)), which was documented in
> {% cite Fisher2015 %}.
> And this technical note was first published as an appendix to [that paper](https://pdfs.semanticscholar.org/396c/b9f172cb681421ed78325a2237bfb428eece.pdf).
> The [FATES documentation](https://fates-docs.readthedocs.io/en/latest/index.html) will provide some more insight on FATES too.
>
{:  .comment}

## Step 1: Get CLM-FATES input data

Preparing CLM-FATES input data is out of scope for this tutorial. We assume the input data tarball contains the following folders:

```
atm   cpl   lnd   share
```

Each sub-folder will then contain all the necessary inputs for running your CLM-FATES case.
For the purpose of this tutorial, input data for a single point location ALP1 (61.0243N, 8.12343E) has been prepared and is ready to use.

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial. If you are not inspired, you can name it *fates*.
>    {% include snippets/create_new_history.md %}
> 2. Import the [input data](https://doi.org/10.5281/zenodo.4108341) and the [restart dataset from Zenodo](https://doi.org/10.5281/zenodo.4126404)
>    or from the shared data library
>
>    ```
>    https://zenodo.org/record/4108341/files/inputdata_version2.0.0_ALP1.tar
>    https://zenodo.org/record/4126404/files/CTSM_FATES-EMERALD_version2.0.0_ALP1_restart_2300-01-01.tar
>    ```
>
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
> 3. Check the datatype (for both files) is **tar**
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
>
> 4. Rename datasets as
>
>    As `https://zenodo.org/record/4108341/files/inputdata_version2.0.0_ALP1.tar` and `https://zenodo.org/record/4126404/files/CTSM_FATES-EMERALD_version2.0.0_ALP1_restart_2300-01-01.tar`
>    are not beautiful names and can give errors for some tools,
>    it is a good practice to change the dataset names by something more meaningful. For example by removing
>    `https://zenodo.org/record/4108341/files/` and `https://zenodo.org/record/4126404/files/` to obtain `inputdata_version2.0.0_ALP1.tar`
>    and `CTSM_FATES-EMERALD_version2.0.0_ALP1_restart_2300-01-01.tar`, respectively.
>
>    {% include snippets/rename_dataset.md %}
>
> 5. Add a tag to the dataset corresponding to `fates`
>
>    {% include snippets/add_tag.md %}
>
{: .hands_on}

# Step 2: Setting up a CLM-FATES simulation

We will be using the CTSM/FATES-EMERALD Galaxy tool.

> ### {% icon comment %} Tip: Finding your tool
>
> Different Galaxy servers may have tools available under different sections, therefore it is often useful to use the **search bar** at the top of the tool panel to find your tool.
>
> Additionally different servers may have multiple, similarly named tools which accomplish similar functions. When following tutorials, you should use precisely the tools that they describe. For real analyses, however, you will need to search among the various options to find the one that works for you.
>
{: .comment}

> ### {% icon hands_on %} Hands-on: Creating a new CTSM/FATES-EMERALD case
>
> 1. {% tool [CTSM/FATES-EMERALD](toolshed.g2.bx.psu.edu/repos/climate/ctsm_fates/ctsm_fates/2.0.1) %} with the following parameters:
>    - {% icon param-file %} *"inputdata for running FATES EMERALD"*: select the **inputdata_version2.0.0_ALP1.tar** file from your history
>    - *Name of your case*: ALP1_exp
>    - Expand *'Customize the model run period'* to change the default values:
>        - **Determines the model run initialization type**: select **hybrid**
>        - **Reference case for hybrid or branch runs**: ALP1_refcase
>        - **Reference date for hybrid or branch runs (yyyy-mm-dd)**: 2300-01-01
>        - **Run start date (yyyy-mm-dd). Only used for startup or hybrid runs**: 0001-01-01
>        - **Restart for running FATES EMERALD**: CTSM_FATES-EMERALD_version2.0.0_ALP1_restart_2300-01-01.tar
>        - **Provides a numerical count for STOP_OPTION**: 5
>        - **Sets the run length along with STOP_N and STOP_DATE**: nyears
>    - Click **Execute**
>
>    > ### {% icon comment %} Tip: search for the tool
>    >
>    > Use the **tools search box** at the top of the tool panel to find **Remove beginning** {% icon tool %}.
>    {: .tip}
>
>    > ## {% icon comment %} startup versus hybrid
>    >
>    >  When using **startup**, the FATES model will start from some arbitrary baseline state that is not linked to any previous run.
>    > Startup runs are typically initialized using a start date of 0001-01-01 except if you change it (start date option).
>    > For any scientific study, starting from an arbitraty baseline state implies you would need to run the model for a long period (between 100 and 200 years)
>    > before being able to use the model outputs. For this reason, we usually make a first simulation (spin-up) in **startup** mode and reuse this case as a baseline
>    > for our scientific study. We then use **hybrid** type and give additional inputs (restart files) to our simulation case. It is then important to specify the dates
>    > of your restart files. This is what we do in this tutorial.
>    >
>    {: .comment}
>
> 2. Check that the datatype of your outputs (history file) is **netcdf**
>
>    All the history files contain gridded data values written at specified times during the model run.
>    Depending on the length of your simulation, you may have one or more history files that you can recognize from their names:
>    `ALP1_exp.clm2.h0.yyyy-mm-dd-sssss.nc` (for non-monthly history files).
>    Datatypes are, by default, automatically guessed. Here, as the prefix is `.nc`, the format is not always recognized as `netcdf` files.
>    To cope with that, one can change the datatype manually, as shown below.
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
>
> 3. **Rename** {% icon galaxy-pencil %} the output dataset (history file) to `ALP1_exp.nc`
>
>    Our FATES model has run for 5 years only, so we get a single output file. As previously, we recommend
>    to rename all netCDF files so that they do not contain any special characters or dots (except for the file extension) or slashes. Some tools, in
>    particular Panoply, won't be able to recognize your file if not named properly.
>
>    {% include snippets/rename_dataset.md %}
>
> 4. Getting metadata information for CLM-FATES netCDF outputs
>
>    {% tool [NetCDF xarray Metadata Info](toolshed.g2.bx.psu.edu/repos/ecology/xarray_metadata_info/xarray_metadata_info/0.15.1) %}  with the following parameters:
>    - **Netcdf file**: ALP1_exp.nc
>    - Click **Execute**
>
>    Inspect the generated output files and identify which variables would provide you some insights about canopy transpiration.
>
>    > ### {% icon question %} Questions
>    >
>    > 1. What are the short names of the relevant variables? Which one will you pick if you want a result in **mm/s**?
>    > 2. What are the dimensions of these variables?
>    >
>    > > ### {% icon solution %} Solution
>    > > 1. **FCTR** is the canopy transpiration in W/m^2 and **QVEGT** is in mm/s. Therefore, we would select the latter.
>    > > 2. These variables are stored as a function of time and lndgrid and since we have only one grid cell, lngrid=1, hence the time series.
>    > {: .solution}
>    {: .question}
>
{: .hands_on}

# Step 3: Quick visualization with Panoply

## Opening up Panoply

> ### {% icon hands_on %} Hands-on: Launch Panoply
>
>  Panoply is available as a Galaxy interactive environment and may not be available on all Galaxy servers.
>
> > ### {% icon tip %} Tip: Launch Panoply in Galaxy
> > Currently Panoply in Galaxy is available on useGalaxy.eu instance, on the "Interactive tools" tool panel section or,
> > as all interactive tools, from the dedicated useGalaxy.eu subdomain: [Live.useGalaxy.eu](https://live.usegalaxy.eu).
> > You may have to login again to [Live.usrGalaxy.eu](https://live.usegalaxy.eu) (use the same username and password than on other useGalaxy.eu subdomains)
> > and switch to the correct history.
> >
> > 1. Open the Panoply tool {% icon tool %} by clicking [here](https://live.usegalaxy.eu/?tool_id=interactive_tool_panoply){:target="_blank"}
> > 2. Check **ALP1_exp.nc** dataset selected in the netcdf input field
> > 3. Click Execute
> > 4. The tool will start running and will stay running permanently
> > 5. Click on the "User" menu at the top and go to "Active Interactive Tools" and locate the Panoply instance you started.
> > 6. Click on your Panoply instance
> >    ![Panoply dataset selection](../../images/select_dataset.png "Select dataset")
> > 7. Click on **ALP1_exp.nc** dataset
> {: .tip}
{: .hands_on}

## Inspect metadata

> ### {% icon hands_on %} Hands-on: Inspect dataset
>
> 1. Inspect dataset content
>
>    Here you can look at the dataset (`ALP1_exp.nc`) and related variables (FSDS, FSA, AREA_TREE, BIOMASS_CANOPY, etc.)
>
>    > ### {% icon question %} Question
>    >
>    > 1. What is the long name of **MORTALITY**?
>    > 2. What is its physical unit?
>    >
>    > > ### {% icon solution %} Solution
>    > > 1. Rate of total mortality per PFT
>    > > 2. indiv/ha/yr
>    > {: .solution}
>    {: .question}
>
>
> 2. Plot the total carbon in live plant leaves (**LEAFC**)
>
>    Cutomize your plot and save it as **png** file in the output folder. Remember that
>    if you do not save in the output folder, your plot will get lost.
>    > ### {% icon question %} Question
>    > 1. Can you observe any pattern? Does it make any sense?
>    >
>    > > ### {% icon solution %} Solution
>    > > 1. We can clearly see a seasonal cycle.
>    > > ![Panoply LEAFC timeserie](../../images/panoply_LEAFC_ALP1_exp.png "LEAFC")
>    > {: .solution}
>    {: .question}
>
>
> 3. Plot the rate of total mortality per PFT (MORTALITY)
>
>    Select a 2D plot with time as x-axis and colored by the rate of total mortality per PFT.
>    Make sure to adjust the y-axis and save your plots in the output folder (as png file).
>    > ### {% icon question %} Question
>    > 1. Can you observe any pattern? Does it make any sense?
>    >
>    > > ### {% icon solution %} Solution
>    > > 1. We can clearly see a seasonal cycle.
>    > > ![Panoply MORTALITY per PFT](../../images/panoply_MORTALITY_ALP1_exp.png "total mortality per PFT")
>    > {: .solution}
>    {: .question}
>
>    > ## {% icon comment %} Quit Panoply properly to save your plots!
>    >
>    > To make sure all your plots stored in **outputs** folder  get exported to Galaxy, you need to quit panoply: **File** --> **Quit Panoply**.
>    {: .comment}
{: .hands_on}

# Step 4: Using Galaxy tools for analysing your CLM-FATES simulation

Panoply is a great tool for exploring the results of your simulations but what we would like is to automate the generation of the plots
so that we can reuse it for any simulations.

> ### {% icon hands_on %} Hands-on: Select and plot **LEAFC**
>
> 1. Select the total carbon in live plant leaves (**LEAFC**)
>    {% tool [NetCDF xarray Selection](toolshed.g2.bx.psu.edu/repos/ecology/xarray_select/xarray_select/0.15.1) %} with the following parameters:
>    - **Input netcdf file**: ALP1_exp.nc
>    - **Tabular of variables**: Metadata infos from ALP1_exp.nc
>    - **Choose the variable to extract**: LEAFC
>    - Click **Execute**
>
> 2. Rename Dataset to **NetCDF xarray Selection on ALP1_exp.nc**
>
>    {% include snippets/rename_dataset.md %}
>
> 3. Clean date column for plotting
>    {% tool [Replace parts of text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_find_and_replace/1.1.3) %} with the following parameters:
>    - **File to process**: NetCDF xarray Selection on ALP1_exp.nc
>    - **Find pattern**:  00:00:00
>    - **Find-Pattern is a regular expression**: Select *No*
>    - **Replace all occurences of the pattern**: Select *Yes*
>    - **Case-Insensitive search**: Select *No*
>    - **Find whole-words**: Select *Yes*
>    - **Ignore first line**: Select *Yes*
>    - **Find and Replace text in**: Select *entire line*
>    - Click **Execute**
>
> 4. Rename Dataset to **LEAFC_clean.tabular**
>
>    {% include snippets/rename_dataset.md %}
>
> 5. Plot the total carbon in live plant leaves (**LEAFC**)
>
> To make a plot, you can use **Scatterplot w ggplot2**  {% icon tool %} with the following parameters:
>    - *"Input in tabular format"*: `LEAFC_clean.tabular`
>    - *"Column to plot on x-axis"*: 1
>    - *"Column to plot on y-axis"*: 4
>    - *"Plot title"*: Total carbon in live plant leaves
>    - *"Label for x axis"*: Time
>    - *"Label for y axis"*: LEAFC (kgC ha-1)
>    - In `Advanced Options` change `Type of plot` to **Points and Lines**.
>    - And finally in `Output options` set `width of output` to **19.0 and `height of output`* to 5.0*.
>
> 6. Click on **Execute**.
>
> 7. **View** {% icon galaxy-eye%} the resulting plot:
>
>    ![LEAFC](../../images/LEAFC_ALP1_exp_ggplot.png)
>
{: .hands_on}

# Step 5: Convert your analysis history into a Galaxy workflow

> ### {% icon hands_on %} Hands-on: Extract workflow
>
> 1. Go to the **History Options menu**  {% icon galaxy-gear %} menu
>    - Select the **Extract Workflow** option.
>    - Remove any unwanted steps
>
> 2. **Rename** the workflow to something descriptive
>    - For example: `CLM-FATES_ ALP1 simulation (5 years)`.
>    - If there are any steps that shouldn't be included in the workflow, you can **uncheck** them.
>
> 3. Click "Create Workflow"
>    - Click on "edit" and check your workflow
>    - Check all the steps
>
{: .hands_on}

# Step 6: Change your CLM-FATES case and rerun your workflow

We would like to run a CLM-FATES case where the atmospheric Carbon Dioxyde Concentration (CO2) is increased by a factor of 4.

> ### {% icon hands_on %} Hands-on: Compare the two simulations
>
>    Using the results from your two CLM-FATES simulations and the generated plots, assess the impact
>    of an increase in the atmosperhic CO2 on the outputs of the model.
>    1. Edit your workflow and customize it to run your new CO2 experiment. For this you would need to
>       add an extra step to extract the first history file from the history collection and generate the
>       corresponding plot. The final workflow would be similar to the one shown below:
>
>  ![FATES workflow](../../images/fates_workflow.png "FATES workflow")
>
>    > ### {% icon question %} Question
>    > 1. Is the model response to this significant increase of atmospheric CO2 what you expected?
>    >   Justify your answer.
>    > 2. Is the current workflow (in particular the variables selected for the plots) the best choice?
>    >   What changes/additions would you recommend?
>    >
>    > > ### {% icon solution %} Solution
>    > > 1. Running 5 years is already sufficient to highlight significant changes.
>    > > ![LEAFC 4xCO2](../../images/LEAFC_ALP1_4CO2_ggplot.png)
>    > > 2. Many suggestions can be given here. One simple addition can be the generation of plots where
>    > > both simulations are represented on the same plot.
>    > {: .solution}
>    {: .question}
>
{: .hands_on}

# Share your work

One of the most important features of Galaxy comes at the end of an analysis. When you have published striking findings, it is important that other researchers are able to reproduce your in-silico experiment. Galaxy enables users to easily share their workflows and histories with others.

To share a history, click on the {% icon galaxy-gear %} icon in the history panel and select `Share or Publish`. On this page you can do 3 things:

1. **Make History Accessible via Link**. This generates a link that you can give out to others. Anybody with this link will be able to view your history.
2. **Make History Accessible and Publish**. This will not only create a link, but will also publish your history. This means your history will be listed under `Shared Data → Histories` in the top menu.
3. **Share with a user**. This will share the history only with specific users on the Galaxy instance.

> ### {% icon comment %} Permissions
> Different servers have different default permission settings. Some servers create all of your datasets completely private to you, while others make them accessible if you know the secret ID.
>
> Be sure to select **Also make all objects within the History accessible** whenever you make a history accessible via link, otherwise whomever you send your link to might not be able to see your history.
{: .comment}

> ### {% icon hands_on %} Hands-on: Share history
>
> 1. Share your history with your neighbour (ask for his/her galaxy username).
> 2. Find the history shared by your neighbour. Histories shared with specific users can be accessed by those users under their top masthead "User" menu under `Histories shared with me`.
{: .hands_on}


> ### {% icon comment %} Publish your history to https://workflowhub.eu/
> One step further is to share your workflow on [https://workflowhub.eu](https://workflowhub.eu) where it
> will be stored in a Galaxy workflow format as well as in [Common Workflow Language](https://www.commonwl.org/).
> It provides standardised workflow identifiers and descriptions needed for workflow discovery, reuse, preservation, interoperability and monitoring and metadata harvesting using standard protocols.
> Please note that [https://workflowhub.eu](https://workflowhub.eu) is still under active development.
{:  .comment}

# Conclusion

{:.no_toc}

We have learnt to run single-point simulations with FATES-CLM and generate workflows for multi-site scenarios.
