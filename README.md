## OVERVIEW

<details>
<summary>1 Summary </summary>
  <ul>
    
  - This repository contains the source code for long- and short-read analyses as published on:
    - Prywes N, Philips NR, Oltrogge LM, Lindner S, Candace Tsai YC, de Pins B, Cowan AE, Taylor-Kearney LJ, Chang HA, Hall LN, Bellieny-Rabelo D, Nisonoff HM, Weissman RF, Flamholz AI, Ding D, Bhatt AY, Shih PM, Mueller-Cajar O, Milo R, Savage DF. A map of the rubisco biochemical landscape. bioRxiv [Preprint]. 2024 Apr 11:2023.09.27.559826. doi: 10.1101/2023.09.27.559826. PMID: 38645011; PMCID: PMC11030240.
  
  </ul>
  
</details>

## GETTING STARTED

<details>
<summary>2 Dependencies</summary>
<ul>

<details>
<summary>2.1 Anaconda </summary>
<ul>

 - Install Miniconda:
 - Download the installer at: https://docs.conda.io/projects/miniconda/en/latest/
 
   ```
   bash Miniconda3-latest-<your-OS>.sh
   ```
  - Set up and update conda: 
    ```
    conda update --all
    conda config --set channel_priority strict
    ```
</ul>
</details>

<details>
<summary>2.2 Snakemake </summary>
<ul>

- Snakemake can be installed directly via Anaconda:

  ```
  conda install -n base -c conda-forge mamba
  ```
</ul>
</details>

</ul>
</details>


<details>
<summary>3 Snakedali Pipeline Setup</summary>
<ul>

<details>
<summary>3.1 Config Files </summary>

  <ul> 

 - Each run can be customized based on the `configuration files`: 
    - `config/pacbio_read_processing.yaml`
    - `config/nextseq_read_processing.yaml`
    
 - From the configuration file users are expected to set up:
   - In-/Output paths for the run
   - Sample names
   - Flanking sequence
  
  </ul>
</details> 

<details>
<summary>3.2 Pipeline Execution </summary>
  <ul>
  
  - This pipeline was designed to work with the `SLURM` job scheduler with the following paramaters:

  - Pacbio Reads
```
snakemake --snakefile pacbio_read_processing.smk -j 5 --cluster "sbatch -t {cluster.time} -n {cluster.cores}" --cluster-config config/cluster.yaml --latency-wait 120 --rerun-incomplete --use-conda 
```
  - NextSeq Reads

```
snakemake --snakefile nextseq_read_processing.smk -j 5 --cluster "sbatch -t {cluster.time} -n {cluster.cores}" --cluster-config config/cluster.yaml --latency-wait 120 --rerun-incomplete --use-conda
```

  - Make sure to adjust the parameters above according to the house rules of your HPC.
</ul>
</details>

</ul>
</details>
