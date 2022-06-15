import scanpy as sc
import sfaira
from typing import List
import argparse


def curate_adata(h5ad_path: str,
                 yaml_path: str,
                 outdir: str = None,
                 custom_fields: List[str] = ['sorting_protocol'],
                 save: bool = True
                 ):
    '''
    Save curated anndata object based on annotation stored in sfaira dataloader .yaml

    Params:
    ------
    - h5ad_path: path to h5ad object
    - yaml_path: path to .yaml file for curation
    - outdir: path to output directory to store curated anndata (default: same directory as input, only used if save=True) 
    - custom_fields: List of custom YAML fields to add (extra from fields in `sfaira.consts.AdataIdsSfaira`)
    - save: boolean indicating whether the curated anndata should be saved in h5ad format or if the function should return the anndata object

    Returns:
    -------
    None, a new h5ad file is saved 
    '''
    # Load anndata object
    print(f"Loading {h5ad_path}...")
    adata_sfaira = sc.read_h5ad(h5ad_path)

    # Make sfaira data object
    d = sfaira.data.interactive.DatasetInteractive(adata_sfaira)

    # Add attributes for custom slots in yaml
    for k in custom_fields:
        setattr(d, k, None)
        setattr(d, f'{k}_obs_key', None)

    # Load annotations in yaml
    print(f"Loading metadata from {yaml_path}...")
    yaml_vals = sfaira.data.utils.read_yaml(fn=yaml_path)
    # Set organism first as this is required to disambiguate valid entries for other meta data.
    k = "organism"
    v = yaml_vals["attr"]["organism"]
    setattr(d, k, v)
    for k, v in yaml_vals["attr"].items():
        if v is not None and k not in ["organism", "sample_fns", "dataset_index"]:
            if isinstance(v, dict):  # v is a dictionary over file-wise meta-data items
                assert d.sample_fn in v.keys(
                ), f"did not find key {d.sample_fn} in yamls keys for {k}"
                v = v[d.sample_fn]
            # Catches spelling errors in meta data definition (yaml keys).
            if not hasattr(d, k) and not hasattr(d, "_" + k):
                raise ValueError(f"Tried setting unavailable property {k}.")
            try:
                setattr(d, k, v)
            except AttributeError as e:
                raise ValueError(
                    f"An error occured when setting {k} as {v}: {e}")

    d.streamline_metadata(keep_orginal_obs=True, clean_obs=False, clean_obs_names=False)

    # Add extra fields to uns
    for k in custom_fields:
        if hasattr(d, k) and getattr(d, k) is not None:
            val = getattr(d, k)
        elif hasattr(d, f"{k}_obs_key") and getattr(d, f"{k}_obs_key") is not None:
            val = d.adata.obs[getattr(d, f"{k}_obs_key")].unique().astype(
                str).tolist()
        else:
            val = 'unknown'
        d.adata.uns[k] = val

    # Add extra fields to obs
    for k in custom_fields:
        if hasattr(d, f"{k}_obs_key") and getattr(d, f"{k}_obs_key") is not None:
            val = d.adata.obs[getattr(d, f"{k}_obs_key")].values.tolist()
            d.adata.obs[k] = val
        else:
            val = 'unknown'

    # d.streamline_var()
    if save:
        if outdir is None:
            out_h5ad_path = h5ad_path.split(".h5ad")[0] + '.curated.h5ad'
        else:
            out_h5ad_path = outdir + \
                h5ad_path.split(".h5ad")[0].split('/')[-1] + '.curated.h5ad'
        print(f"Saving curated data to {out_h5ad_path}...")
        adata_sfaira.write_h5ad(out_h5ad_path)
    else:
        return(adata_sfaira)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "h5ad_path", help="path to h5ad object containing dataset")
    parser.add_argument(
        "yaml_path", help="path to yaml object containing metadata")
    parser.add_argument("--outdir",
                        default=None,
                        help="folder to save curated h5ad (uses same as input by default)")
    args = parser.parse_args()
    curate_adata(args.h5ad_path, args.yaml_path, outdir=args.outdir)
