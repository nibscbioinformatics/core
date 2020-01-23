# Analysing the CellRanger data with Serurat


This documents explains, with examples, how to call each R script to analyze scRNA-seq data, using demultiplexed files from cellranger.
The files 'aggregate_example.html', 'analyse_unfiltered_example.html' and 'analyse_filtered_example.html' are examples of each of the following scripts.

## Aggregate the samples

	Aggregate samples:
		```
    Rscript -e "rmarkdown::render('aggregate.Rmd', params = list(sample_paths = '/path/to/sample_1,/path/to/sample_2,/path/to/sample_3', sample_names = 'sample_name_1,sample_name_2,sample_name_3', output_path = '/path/to/output'))"
    ```

		Example:
		```
    Rscript -e "rmarkdown::render('aggregate.Rmd', params = list(sample_paths = '/home/AD/praposo/software/cell_ranger/cellranger-3.0.2/MSC/MSC_components/bone_marrow/outs/filtered_feature_bc_matrix/,/home/AD/praposo/software/cell_ranger/cellranger-3.0.2/MSC/MSC_components/bone_marrow_mb/outs/filtered_feature_bc_matrix/,/home/AD/praposo/software/cell_ranger/cellranger-3.0.2/MSC/MSC_components/adipose/outs/filtered_feature_bc_matrix/', sample_names = 'BM,BM2,Adipose', output_path = '/home/AD/praposo/WGS/'))"
		```

		Note: ```/path/to/sample/``` is the cellranger output folder ```/outs/filtered_feature_bc_matrix/```

## Analyse unfiltered data

	Analysis (unfiltered):
		```
    Rscript -e "rmarkdown::render('analyse_unfiltered.Rmd', params = list(input_path = '/path/to/input'))"
		```

		Example:
		```
    Rscript -e "rmarkdown::render('analyse_filtered.Rmd', params = list(input_path = '/home/AD/praposo/WGS/scRNA_nextflow/'))"
		```

		Note: ```/path/to/input``` is folder where ```'aggregated_object'``` is found

## Analyse the filtered data

	Analysis (filtered):
		```
    Rscript -e "rmarkdown::render('analyse_filtered.Rmd', params = list(input_path = '/path/to/input'))"
		```

		Example:
		```
    Rscript -e "rmarkdown::render('analyse_filtered.Rmd', params = list(input_path = '/home/AD/praposo/WGS/scRNA_nextflow/'))"
		```
    
		Note: /path/to/input is folder where 'aggregated_object' is found
