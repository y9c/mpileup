# mpileup

## Install

```bash
cargo install mpileup
```

## Usage

count number of base and indel

```bash
mpileup --target test/region.bed --reference test/reference.fa --input test/sample1.bam test/sample2.bam -g -d 10
```

In this examples,

- pileup up reads within region in `region.bed`
- use `reference.fa` as reference
- accept multiple input bam files: `sample1.bam`, `sample2.bam` ...
- report indel with argument `-g`
- set depth cutoff as 10 by `-d 10`
