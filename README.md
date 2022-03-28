# mpileup

## Install

```bash
cargo install mpileup
```

## Usage

### pileup number of base and indel

```bash
mpileup base --target test/region.bed --reference test/reference.fa --input test/sample1.bam test/sample2.bam -g -d 10
```

In this subcommand (example),

- pileup up reads within region in `region.bed`
- use `reference.fa` as reference
- accept **multiple** input bam files: `sample1.bam`, `sample2.bam` ...
- report indel with argument `-g`
- set depth cutoff as 10 by `-d 10`

### count number of reads

```bash
mpileup count --target test/region.bed --reference test/reference.fa --input test/sample1.bam
```

In this subcommand (example),

- only support one file a time

## TODO

- support sam flag filtering
- support trimming read ends
- remove duplicate by UMI on qname
- group the result by customized tag (such as cell barcode)
