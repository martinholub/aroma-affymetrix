---
title: Config Documentation
output: html_document
---

# Config for `aroma` big experiments processing

The aroma pipeline is completely controlled by the config. This is in contrast to the classic $\mu$Array where many arguments are passed on command line. The pipeline was successfully tested on two rat and one human experiment with around 8000 samples each.

The pipeline is executed with:
```bash
Rscript --no-save --verbose ma_bigexp_norm.R path/to/config_aroma.yml
# Rscript --no-save --verbose ~/git/git-dev/QC_Normalization_V4/maQCN_pipeline/ma_bigexp_norm.R ~/git/git-dev/QC_Normalization_V4/tests/config_aroma.yml
# Rscript --no-save --verbose ~/QC_Normalization_V4/maQCN_pipeline/ma_bigexp_norm.R /rnaseq/results/NEB_P1/RN/RN-00187aroma/config_aroma.yml
```
If the config is not passed on command line, the one in `configs` folder of the $\mu$Array pipeline is used. The config looks something like this:

```yml
# Aroma Affymetrix Big experiment processing config
default:
  raw_path: "/biodata/ma-raw/RN-00186" # folder with CEL files (can be CEL.gz)
  exp_path: "/rnaseq/results/NEB_P1/RN/RN-00186" # parent folder of processing and normalization results
  # chipdb can be also name of installed annotation package or a path to ASCII CDF file downloaded from BA website, e.g.
  #chipdb: "~/tmp/big_exp/test/annotationData/chipTypes/rat2302rnensg/rat2302rnensg.cdf"
  chipdb: "rat2302cdf"
  exp_name: # if exp_name is empty, it will be pulled from raw_path
  rmaTrim: 0.05 # dont change
  rmaTv: 1000 # dont change
```

## raw_path:
Describes location where the raw data are stored as gzipped `CEL` files.

## exp_path
Describes location where the aroma folder structure will be created. This location holds all intermediate files as well as results.

Note that the intermediate files are used on subsequent reprocessing if available. This potentially results in skipping most of the processing.

## chipdb
This parameter is either name of installed annotation package or path to an ASCII CDF file (e.g. downloaded from BrainArray website). Also binary CDF files are accepted. In any case, the CDF used for the processing is a binary CDF file so variable number of steps must be taken to obtain it.

## exp_name
The NBR epxeriment name. Can be left empty, then it will be pulled from `raw_path`. Must be specified if experiment has different name in `exp_path`.

## rmaTrim, rmaTv
**Do not change.** These values must correspond to the values used my classic $\mu$Array
