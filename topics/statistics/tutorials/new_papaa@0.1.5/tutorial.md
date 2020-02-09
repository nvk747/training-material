---
layout: tutorial_hands_on

title: Pancancer Aberrant Pathway Activity Analysis
zenodo_link: https://zenodo.org/record/3632117#.Xjh8MpNKi-U
questions:
- how to predict aberrant pathway activity in The cancer genome Atlas(TCGA) using
  Machine learning approaches?
objectives:
- Learn to predict aberrant pathway activity using RNAseq data, mutational status
  and copy number variation data from TCGA.
- Apply logistic regression based machine learning algorithms on TCGA data.
time_estimation: ''
key_points:
- Identify the transcriptomic signature of tumor and potential biomarkers that are
  useful in improving personalized medicine
- They will appear at the end of the tutorial
contributors:
- nvk747
- blankenberg

---


# Introduction
{:.no_toc}

<!-- This is a comment. -->

General introduction about the topic and then an introduction of the
tutorial (the questions and the objectives). It is nice also to have a
scheme to sum up the pipeline used during the tutorial. The idea is to
give to trainees insight into the content of the tutorial and the (theoretical
and technical) key concepts they will learn.

You may want to cite some publications; this can be done by adding citations to the
bibliography file (`tutorial.bib` file next to your `tutorial.md` file). These citations
must be in bibtex format. If you have the DOI for the paper you wish to cite, you can
get the corresponding bibtex entry using [doi2bib.org](https://doi2bib.org).

With the example you will find in the `tutorial.bib` file, you can add a citation to
this article here in your tutorial like this:
{% raw %} `{% cite Batut2018 %}`{% endraw %}.
This will be rendered like this: {% cite Batut2018 %}, and links to a
[bibliography section](#bibliography) which will automatically be created at the end of the
tutorial.


**Please follow our
[tutorial to learn how to fill the Markdown]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-content/tutorial.html)**

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Title for your first section

Give some background about what the trainees will be doing in the section.
Remember that many people reading your materials will likely be novices,
so make sure to explain all the relevant concepts.

## Title for a subsection
Section and subsection titles will be displayed in the tutorial index on the left side of
the page, so try to make them informative and concise!

# Hands-on Sections
Below are a series of hand-on boxes, one for each tool in your workflow file.
Often you may wish to combine several boxes into one or make other adjustments such
as breaking the tutorial into sections, we encourage you to make such changes as you
see fit, this is just a starting point :)

Anywhere you find the word "***TODO***", there is something that needs to be changed
depending on the specifics of your tutorial.

have fun!

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]() or from the shared data library
>
>    ```
>    
>    ```
>    ***TODO***: *Add the files by the ones on Zenodo here (if not added)*
>
>    ***TODO***: *Remove the useless files (if added)*
>
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
> 3. Rename the datasets
> 4. Check that the datatype
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
>
> 5. Add to each database a tag corresponding to ...
>
>    {% include snippets/add_tag.md %}
>
{: .hands_on}

# Title of the section usually corresponding to a big step in the analysis

It comes first a description of the step: some background and some theory.
Some image can be added there to support the theory explanation:

![Alternative text](../../images/image_name "Legend of the image")

The idea is to keep the theory description before quite simple to focus more on the practical part.

***TODO***: *Consider adding a detail box to expand the theory*

> ### {% icon details %} More details about the theory
>
> But to describe more details, it is possible to use the detail boxes which are expandable
>
{: .details}

A big step can have several subsections or sub steps:


## Sub-step with **PAPAA: PanCancer classifier**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **PAPAA: PanCancer classifier** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Filename of features to use in model"*: `output` (Input dataset)
>    - *"Treat X matrix as 'raw'"*: `Yes`
>    - {% icon param-file %} *"Filename mutations"*: `output` (Input dataset)
>    - {% icon param-file %} *"Filename of mutation burden"*: `output` (Input dataset)
>    - {% icon param-file %} *"Filename of sample"*: `output` (Input dataset)
>    - *"Comma separated string of HUGO gene symbols"*: `ERBB2,PIK3CA,KRAS,AKT1`
>    - *"Comma sep string of TCGA disease acronyms. If no arguments are passed, filtering will default to options given in --filter_count and --filter_prop."*: `BLCA,BRCA,CESC,COAD,ESCA,LUAD,LUSC,OV,PRAD,READ,STAD,UCEC,UCS`
>    - *"Decision to drop input genes from X matrix"*: `Yes`
>    - *"Supplement Y matrix with copy number events"*: `Yes`
>    - *"the alphas for parameter sweep"*: `0.1,0.13,0.15,0.18,0.2,0.3,0.4,0.6,0.7`
>    - *"the l1 ratios for parameter sweep"*: `0.1,0.125,0.15,0.2,0.25,0.3,0.35`
>    - *"alternative genes to test performance"*: `PTEN,PIK3R1,STK11`
>    - *"The alternative diseases to test performance"*: `BRCA,COAD,ESCA,HNSC,LGG,LUAD,LUSC,PRAD,READ,GBM,UCEC,UCS`
>    - *"Remove hypermutated samples"*: `Yes`
>    - *"Keep intermediate ROC values for plotting"*: `Yes`
>    - *"Shuffle the input gene exprs matrix alongside"*: `Yes`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **PAPAA: PanCancer within disease analysis**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **PAPAA: PanCancer within disease analysis** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Filename of features to use in model"*: `output` (Input dataset)
>    - *"Treat X matrix as 'raw'"*: `Yes`
>    - {% icon param-file %} *"Filename mutations"*: `output` (Input dataset)
>    - {% icon param-file %} *"Filename of mutation burden"*: `output` (Input dataset)
>    - {% icon param-file %} *"Filename of sample"*: `output` (Input dataset)
>    - *"Comma separated string of HUGO gene symbols"*: `ERBB2,PIK3CA,KRAS,AKT1`
>    - *"Comma sep string of TCGA disease acronyms. If no arguments are passed, filtering will default to options given in --filter_count and --filter_prop."*: `BLCA,BRCA,CESC,COAD,ESCA,LUAD,LUSC,OV,PRAD,READ,STAD,UCEC,UCS`
>    - {% icon param-file %} *"File with Copy number loss"*: `output` (Input dataset)
>    - {% icon param-file %} *"File with Copy number gain"*: `output` (Input dataset)
>    - {% icon param-file %} *"File with cancer gene classification table"*: `output` (Input dataset)
>    - *"the alphas for parameter sweep"*: `0.1,0.13,0.15,0.18,0.2,0.3,0.4,0.6,0.7`
>    - *"the l1 ratios for parameter sweep"*: `0.1,0.125,0.15,0.2,0.25,0.3,0.35`
>    - *"Remove hypermutated samples"*: `Yes`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **PAPAA: PanCancer apply weights**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **PAPAA: PanCancer apply weights** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Filename of features to use in model"*: `output` (Input dataset)
>    - {% icon param-file %} *"Filename mutations"*: `output` (Input dataset)
>    - {% icon param-file %} *"Filename of mutation burden"*: `output` (Input dataset)
>    - {% icon param-file %} *"Filename of sample"*: `output` (Input dataset)
>    - *"Supplement Y matrix with copy number events"*: `Yes`
>    - {% icon param-file %} *"pancancer classifier summary"*: `classifier_summary` (output of **PAPAA: PanCancer classifier** {% icon tool %})
>    - {% icon param-file %} *"pancancer classifier coefficients"*: `classifier_coefficients` (output of **PAPAA: PanCancer classifier** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **PAPAA: PanCancer external sample status prediction**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **PAPAA: PanCancer external sample status prediction** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Classifier data"*: `classifier_summary` (output of **PAPAA: PanCancer classifier** {% icon tool %})
>    - {% icon param-file %} *"pancancer classifier coefficients"*: `classifier_coefficients` (output of **PAPAA: PanCancer classifier** {% icon tool %})
>    - {% icon param-file %} *"external sample gene expression data"*: `output` (Input dataset)
>    - {% icon param-file %} *"given mutational status"*: `output` (Input dataset)
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **PAPAA: PanCancer compare within models**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **PAPAA: PanCancer compare within models** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"pancancer classifier summary"*: `classifier_summary` (output of **PAPAA: PanCancer classifier** {% icon tool %})
>    - {% icon param-file %} *"pancancer classifier coefficients"*: `classifier_coefficients` (output of **PAPAA: PanCancer classifier** {% icon tool %})
>    - {% icon param-file %} *"pan_within classifier summary"*: `classifier_summary` (output of **PAPAA: PanCancer within disease analysis** {% icon tool %})
>    - {% icon param-file %} *"pan_within classifier coefficients"*: `classifier_coefficients` (output of **PAPAA: PanCancer within disease analysis** {% icon tool %})
>    - *"Would you want to compare given model with alt gene model?"*: `do not do alt gene`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **PAPAA: PanCancer visualize decisions**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **PAPAA: PanCancer visualize decisions** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"pancancer decisions"*: `classifier_decisions` (output of **PAPAA: PanCancer apply weights** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **PAPAA: PanCancer map mutation class**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **PAPAA: PanCancer map mutation class** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"pancancer decisions"*: `classifier_decisions` (output of **PAPAA: PanCancer apply weights** {% icon tool %})
>    - {% icon param-file %} *"string of the genes to extract or genelist file"*: `output` (Input dataset)
>    - *"Supplement Y matrix with copy number events"*: `Yes`
>    - {% icon param-file %} *"Filename of raw mut MAF"*: `output` (Input dataset)
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **PAPAA: PanCancer alternative genes pathwaymapper**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **PAPAA: PanCancer alternative genes pathwaymapper** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"pancancer decisions"*: `classifier_decisions` (output of **PAPAA: PanCancer apply weights** {% icon tool %})
>    - *"Comma separated string of HUGO gene symbols"*: `ERBB2,PIK3CA,KRAS,AKT1`
>    - {% icon param-file %} *"string of the genes to extract or genelist file"*: `output` (Input dataset)
>    - {% icon param-file %} *"Filename mutations"*: `output` (Input dataset)
>    - {% icon param-file %} *"Filename of sample"*: `output` (Input dataset)
>    - *"Supplement Y matrix with copy number events"*: `Yes`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **PAPAA: PanCancer pathway count heatmaps**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **PAPAA: PanCancer pathway count heatmaps** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"pancancer decisions"*: `classifier_decisions` (output of **PAPAA: PanCancer apply weights** {% icon tool %})
>    - {% icon param-file %} *"pancancer metrics pathwaymapper"*: `pathway_metrics_pathwaymapper` (output of **PAPAA: PanCancer alternative genes pathwaymapper** {% icon tool %})
>    - {% icon param-file %} *"pancancer gene metric ranks"*: `all_gene_metric_ranks` (output of **PAPAA: PanCancer alternative genes pathwaymapper** {% icon tool %})
>    - *"Comma separated string of HUGO gene symbols"*: `ERBB2,PIK3CA,KRAS,AKT1`
>    - {% icon param-file %} *"String of the pathway genes to extract"*: `output` (Input dataset)
>    - {% icon param-file %} *"Filename of features to use in model"*: `output` (Input dataset)
>    - {% icon param-file %} *"Filename mutations"*: `output` (Input dataset)
>    - {% icon param-file %} *"Filename of mutation burden"*: `output` (Input dataset)
>    - {% icon param-file %} *"Filename of sample"*: `output` (Input dataset)
>    - {% icon param-file %} *"File with Copy number loss"*: `output` (Input dataset)
>    - {% icon param-file %} *"File with Copy number gain"*: `output` (Input dataset)
>    - {% icon param-file %} *"File with cancer gene classification table"*: `output` (Input dataset)
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **PAPAA: PanCancer targene summary figures**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **PAPAA: PanCancer targene summary figures** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Classifier data"*: `classifier_summary` (output of **PAPAA: PanCancer classifier** {% icon tool %})
>    - {% icon param-file %} *"pancancer classifier coefficients"*: `classifier_coefficients` (output of **PAPAA: PanCancer classifier** {% icon tool %})
>    - {% icon param-file %} *"summary counts"*: `summary_counts` (output of **PAPAA: PanCancer classifier** {% icon tool %})
>    - {% icon param-file %} *"mutation classification scores"*: `mutation_classification_scores` (output of **PAPAA: PanCancer map mutation class** {% icon tool %})
>    - {% icon param-file %} *"path events per sample"*: `path_events_per_sample` (output of **PAPAA: PanCancer pathway count heatmaps** {% icon tool %})
>    - {% icon param-file %} *"pancancer gene metric ranks"*: `all_gene_metric_ranks` (output of **PAPAA: PanCancer alternative genes pathwaymapper** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **PAPAA: PanCancer targene cell line predictions**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **PAPAA: PanCancer targene cell line predictions** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Classifier data"*: `classifier_summary` (output of **PAPAA: PanCancer classifier** {% icon tool %})
>    - {% icon param-file %} *"pancancer classifier coefficients"*: `classifier_coefficients` (output of **PAPAA: PanCancer classifier** {% icon tool %})
>    - {% icon param-file %} *"nucleotide mutation scores"*: `nucleotide acid mutation scores` (output of **PAPAA: PanCancer targene summary figures** {% icon tool %})
>    - {% icon param-file %} *"amino acid mutation scores"*: `amino acid mutation scores` (output of **PAPAA: PanCancer targene summary figures** {% icon tool %})
>    - *"Comma separated string of HUGO targene symbols"*: `ERBB2_MUT,PIK3CA_MUT,KRAS_MUT,AKT1_MUT`
>    - {% icon param-file %} *"string of the genes to extract or genelist file"*: `output` (Input dataset)
>    - {% icon param-file %} *"Filename ccle rnaseq data"*: `output` (Input dataset)
>    - {% icon param-file %} *"Filename ccle mutational data"*: `output` (Input dataset)
>    - {% icon param-file %} *"Filename ccle variant data"*: `output` (Input dataset)
>    - {% icon param-file %} *"Filename gdsc rnaseq data"*: `output` (Input dataset)
>    - {% icon param-file %} *"Filename gdsc mutational data"*: `output` (Input dataset)
>    - {% icon param-file %} *"Filename for gdsc1 pharmacological data file"*: `output` (Input dataset)
>    - {% icon param-file %} *"Filename for gdsc2 pharmacological data file"*: `output` (Input dataset)
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **PAPAA: PanCancer targene pharmacology**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **PAPAA: PanCancer targene pharmacology** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"gdsc1 targene pharmacology predictions"*: `gdsc1 targene pharmacology predictions` (output of **PAPAA: PanCancer targene cell line predictions** {% icon tool %})
>    - {% icon param-file %} *"gdsc2 targene pharmacology predictions"*: `gdsc2 targene pharmacology predictions` (output of **PAPAA: PanCancer targene cell line predictions** {% icon tool %})
>    - {% icon param-file %} *"gdsc1 ccle targene pharmacology predictions"*: `gdsc1 ccle targene pharmacology predictions` (output of **PAPAA: PanCancer targene cell line predictions** {% icon tool %})
>    - {% icon param-file %} *"gdsc2 ccle targene pharmacology predictions"*: `gdsc2 ccle targene pharmacology predictions` (output of **PAPAA: PanCancer targene cell line predictions** {% icon tool %})
>    - {% icon param-file %} *"Filename list of compounds"*: `output` (Input dataset)
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}


## Re-arrange

To create the template, each step of the workflow had its own subsection.

***TODO***: *Re-arrange the generated subsections into sections or other subsections.
Consider merging some hands-on boxes to have a meaningful flow of the analyses*

# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.