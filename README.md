# Cas13a guide design

## Data needed
- reference genome sequence (.fa)
- sequenced viral genomes representing genetic diversity (.fa)
- Cas13a repeat sequence

## Need to install:
- [bowtie](http://bowtie-bio.sourceforge.net/index.shtml)
- [RNAfold](https://www.tbi.univie.ac.at/RNA/RNAfold.1.html)
- [R](https://www.r-project.org/), with libraries: optparse, here, Biostrings

## Workflow (using `scripts`)

1. Break reference genome into windows of 20nt: `break_genome.R`
- need reference genome .fa file

2. Score window sensitivity to SARS-CoV-2: `score_SARS-CoV-2_sensitivity.R`
- need .fa file of sequenced viral genomes representing genetic diversity

3. Score specificity against other human coronaviruses: `score_human_CoV_specificity.R`
- need .fa file of offtarget coronavirus genomes

4. Calculate crRNA folding structures (RNAfold): `score_RNAfold_crRNAs.R`
- TODO: specify Cas13a repeat sequence in script arguments

5. Calculate base-pairing propensity of viral target window: `score_RNAfold_target.R`

6. Align against hg38 transcriptome: `align_bowtie.R`

7. Align against bosTau9 transcriptome: `align_bowtie.R`

8. Compile scores into aggregate summary table: `compile_scores.R`

**sample processing script**: `SARS-CoV-2_guidedesign/outputs/cas13a_20nt/process.sh`
