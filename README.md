# Semileptonic VV VBS Selection

## How to run (Interactively)
--------

1. **Make sample list**

`make_sample_list.sh` script takes `eos` path as argument and, looks for directories within and all the root files within each directory will be put into one text file.
For example, 
```bash
mkdir 2018
cd 2018
./make_sample_list.sh /eos/uscms/.../2018_parent_sample_dir
cd -
```
will produce

```
/eos/uscms/.../2018_parent_sample_dir/sample1` -> `2018/sample1.txt
/eos/uscms/.../2018_parent_sample_dir/sample2` -> `2018/sample2.txt
...
```

where `sample1.txt` will contain all the `.root` files within `/eos/uscms/.../2018_parent_sample_dir/sample1`

2. **Run**

`run.sh` script setups `ROOT` and runs the macro `vbs_flat_ntupler.cc`

Arguments:
```
./run.sh    SAMPLE_LIST_TXT    YEAR    EOS_DIR
SAMPLE_LIST_TXT -> path to txt file,
    .cc macro will extract filename produce root file in current directory
YEAR -> sample year
EOS_DIR (Optional) -> copies .root file to the eos directory, if given.
                      useful when running on condor
```

## How to run on HTCondor at lpc
--------

1. Make the files list as decribed before.

2. Change the parameters in `condor_submit.jdl` (mainly in starting)

3. Make the `log_dir` as set in `.jdl` file.

4. Submit

```bash
condor_submit condor_submit.jdl
```
