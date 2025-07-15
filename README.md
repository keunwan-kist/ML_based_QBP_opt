# ML‑Guided QBP Optimization Pipeline

A quick‑start workflow that couples Rosetta structure generation with machine‑learning‑guided sequence optimization.

---

## Prerequisites

| Software | Tested Version |
|----------|---------------|
| **Rosetta** (for RosettaScript) | 2020.28.post.dev+102.master.43e678f |
| **Ranger** (R package) | 0.8.0 |
| **PRROC** (R package) | 1.3 |

Install Rosetta as described on the [Rosetta Commons](https://www.rosettacommons.org/) site and the R packages via CRAN:

```r
install.packages(c("ranger", "PRROC"))
```

---

## Example Workflow

All paths below are relative to the repository root.

1. **Move to the working directory**

   ```bash
   cd test_out
   ```

2. **Generate a single mutant structure**

   ```bash
   rosetta_scripts.linuxgccrelease \
		-parser:protocol ../scripts/mut_gen_relax.xml \
		-s ../input_data/1wdna_0001.pdb \
		-ex1 -ex2 -use_input_sc \
		-resfile ../input_data/rand_resfile2442.txt \
		-out:prefix mut. \
		-nstruct 1 \
		-linmem_ig 10
   # Output : mut.1wdna_0001_0001.pdb
   ```

3. **Graft QBP (chain X) onto the mutant**

   ```bash
   cat mut.1wdna_0001_0001.pdb ../input_data/QBP.pdb > complex.mut.1wdna_0001_0001.pdb
   ```

4. **Minimize and score the protein‑ligand complex**

   ```bash
   rosetta_scripts.linuxgccrelease \
		-parser:protocol ../scripts/min_and_scoring_lig_com.xml \
		-s complex.mut.1wdna_0001_0001.pdb \
		-use_input_sc -ex1 -ex2 \
		-extra_res_fa ../input_data/QBP.params \
		-nstruct 1 \
		-linmem_ig 10
   # Output : complex.mut.1wdna_0001_0001_0001.pdb
   # Scores in the PDB file: lig_ddg, lig_sc, score_filter
   ```

5. **Model Training and Saving**

	```bash
	Rscript ../scripts/reg_model_gen.r ../input_data/seq_ddg_values_3nd.feat.mat ddg
	mv ddg.pred_model.ranger ../pred_model/
	# Saved model : ddg.pred_model.ranger
	# Pre‑trained models are available in ./pred_model/
	```

6. **Run GA‑based optimization (population 100, 10 iterations)**

   ```bash
   perl -I ../scripts/ GA_opt_flex_seq.pl \
		-cov_model ../pred_model/ddg.pred_model.ranger \
		-energy_model ../pred_model/energy.pred_model.ranger \
		-N 100 -iter 10 > opt_flex_seq.out 2> log.txt
   ```

---

## Additional Utilities

### Sequence Feature Generation (5 × L)

```bash
Rscript ../scripts/seqfeat_gen.r ../input_data/test.fasta > test.out.mat
```

---

## Folder Structure (minimal)

```
.
├── input_data/       # PDBs, resfiles, QBP params
├── pred_model/       # Pre‑trained ranger models
├── scripts/          # RosettaScripts XML, Perl, R
└── test_out/         # Working directory for example run
```

Feel free to modify paths or parameters to suit your setup.
