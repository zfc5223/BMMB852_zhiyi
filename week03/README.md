# Week 3: Collaboration and Pull Request

## Repository

I reviewed one of the assigned classmate repositories.

- Original repository: `https://github.com/hec5351/Applied-Bioinformatics`
- My fork: `https://github.com/zfc5223/Haley_Applied-Bioinformatics`
- Local clone of my fork: `/home/zhiyi/1data/projects/review_tmp/Haley_Applied-Bioinformatics`
- Review branch: `fixhaley`
- Pull request: `https://github.com/hec5351/Applied-Bioinformatics/pull/1`

I used a fork because I needed my own copy where I could push changes. I cloned my fork locally, made the change on a separate branch, and will contribute the change back to the original repository with a pull request.

## 1. Fork and Clone

I forked Haley's repository on GitHub, then cloned my fork locally:

```bash
mkdir -p /home/zhiyi/1data/projects/review_tmp
cd /home/zhiyi/1data/projects/review_tmp
git clone https://github.com/zfc5223/Haley_Applied-Bioinformatics.git
cd Haley_Applied-Bioinformatics
git checkout -b fixhaley
```

## 2. Security Check

Before running the code, I inspected the repository:

```bash
find . -maxdepth 3 -type f | sort
sed -n '1,240p' Week02Assignment/Makefile
sed -n '1,260p' Week02Assignment/README.md
git status --short --branch
```

I did not find commands that looked dangerous, such as `sudo`, credential access, or executing an uninspected remote script. The `Makefile` downloaded public genome data from NCBI. The managed risk was acceptable after inspection.

## 3. README Evaluation

The Week 2 README described the selected genome, _Apis florea_, and included answers to the assignment questions.

The main problem was that the README and `Makefile` were not fully consistent. The README described downloading FASTA and GFF files, but the original `Makefile` only produced a FASTA file. The README also did not clearly describe the index files needed for IGV.

## 4. Running Instructions and Outcomes

The expected outcome was mostly clear from the README, but the command for reproducing the workflow was not. The README showed:

```bash
make fasta
make gff
```

Those were not targets in the actual `Makefile`. A reviewer would have difficulty reproducing the results directly from the README.

## 5. Reproducibility Check

I first checked what the original `Makefile` would do:

```bash
pixi run -m "$HOME/edu/bioinfo" make -n -C Week02Assignment all
```

The dry run showed that the original workflow only created:

```text
GCF_048593485.1.fasta
```

It did not create:

```text
GCF_048593485.1.gff.gz
GCF_048593485.1.fasta.fai
GCF_048593485.1.gff.gz.tbi
```

I also confirmed that the course pixi environment had the needed tools:

```bash
pixi run -m "$HOME/edu/bioinfo" samtools --version
pixi run -m "$HOME/edu/bioinfo" tabix --help
```

## 6. AI Comparison With My Solution

I asked the AI agent to compare my Week 2 solution with Haley's Week 2 solution.

Both workflows used NCBI as the data source and used a `Makefile` to automate data retrieval. My solution also prepared the FASTA and GFF files for IGV by creating a FASTA index and a tabix-indexed GFF.

The main difference was that Haley's original `Makefile` only generated the FASTA file, while the README described both FASTA and GFF data.

## 7. Which Solution Was Better

The AI agent judged the workflow that creates both data files and index files as better for this assignment. I agree with that evaluation because the assignment asks for both obtaining genomic data and visualizing the genome in IGV.

## 8. Summary of Findings

This was not a problem with the genome choice. _Apis florea_ is a valid genome of interest for the assignment.

The main reproducibility issue was that the README and `Makefile` did not match. A reviewer following the README would expect both FASTA and GFF outputs, but the original `Makefile` only produced FASTA and did not create IGV index files.

## 9. Change Made

I made a small targeted change to the forked repository.

Changed files:

```text
Week02Assignment/Makefile
Week02Assignment/README.md
```

The updated `Makefile` now downloads the _Apis florea_ FASTA and GFF files from NCBI FTP, creates a FASTA index with `samtools faidx`, and creates a tabix-indexed GFF with `bgzip` and `tabix`.

The README was updated so the reproduction command matches the actual workflow:

```bash
cd Week02Assignment/
pixi run -m "$HOME/edu/bioinfo" make all
```

## 10. Commit and Push

The change is ready to commit and push from the fork branch:

```bash
cd /home/zhiyi/1data/projects/review_tmp/Haley_Applied-Bioinformatics
git add Week02Assignment/Makefile Week02Assignment/README.md
git commit -m "Add Week 2 genome indexes"
git push origin fixhaley
```

## 11. Pull Request

After pushing the branch, I will open a pull request from:

```text
zfc5223:fixhaley
```

to:

```text
hec5351:main
```

Pull request URL:

```text
https://github.com/hec5351/Applied-Bioinformatics/pull/1
```

## 12. Author Review

After the pull request is opened, the repository owner can review the change. They may merge it, request changes, or decline it. The fork will keep my branch either way.

## Validation After the Change

I validated the updated workflow with the course pixi environment:

```bash
pixi run -m "$HOME/edu/bioinfo" make -C Week02Assignment all
```

This completed successfully and produced:

```text
GCF_048593485.1.fasta
GCF_048593485.1.fasta.fai
GCF_048593485.1.gff.gz
GCF_048593485.1.gff.gz.tbi
```

I also checked that the index files worked:

```bash
pixi run -m "$HOME/edu/bioinfo" samtools faidx GCF_048593485.1.fasta NW_028256300.1:1-60
pixi run -m "$HOME/edu/bioinfo" tabix GCF_048593485.1.gff.gz NW_028256300.1:1-38558
```

Both commands returned output. After validation, I ran `make clean` so generated genome files would not be included in the pull request.
