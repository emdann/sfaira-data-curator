Curate and harmonize single-cell datasets stored locally in h5ad files. 

This code uses functionality in the [sfaira](https://sfaira.readthedocs.io/en/latest/index.html) library (thanks for sharing the SCOG )

## Why
- Store metadata in a human and machine readable format
- Automate conversion of anndata objects to a common set of columns and slots in .obs, .var, .uns

## Usage instructions

### Install sfaira
```
pip install sfaira
```

### Locate anndata file
```
H5AD_PATH=/path/to/h5ad/file.h5ad
```

### Prepare curation yaml
For each dataset with similar metadata schema, we will collect information about standard metadata in a YAML file. This file uses the annotation format used for sfaira dataloaders (summarised below, see [docs](https://sfaira.readthedocs.io/en/latest/adding_datasets.html#field-descriptions) for detailed information on all the YAML fields).

- **dataset-wise fields:** here we store information at the dataset level, mostly used to make ids
- **layers:** here we indicate where to find raw and processed data matrices in the anndata 
- **dataset_or_feature_wise:** info about features. They can be supplied as `NAME` or as `NAME_var_key`: The former indicates that the entire data set has one value, which is stated in the yaml. The latter, `NAME_var_key`, indicates that there is a column in `adata.var[NAME_var_key]` which contains the annotation per feature for this meta data item.
- **dataset_or_observation_wise:** info about observations (cells). They can be supplied as `NAME` or as `NAME_obs_key`: The former indicates that the entire data set has one value, which is stated in the yaml. The latter, `NAME_obs_key`, indicates that there is a column in `adata.obs[NAME_obs_key]` which contains the annotation per cell for this meta data item. For example, if all the cells come from blood samples, the `organ` annotation should be set in the yaml as:
```
    organ: 'blood'
    organ_obs_key: 
```  
If the dataset collects samples from multiple organs, and the annotation is stored in `adata.obs['organ']`, the `organ` annotation should be set in the yaml as:
```
    organ:
    organ_obs_key: 'organ'
```  


Additional lab specific fields include:
- `sorting_protocol`: which FACS sorting protocol was  

### Save curated dataset
```
curate_adata.py
```
 


