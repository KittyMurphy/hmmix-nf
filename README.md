# Hmmix-nf

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.4-23aa62.svg)](https://www.nextflow.io/)

## Introduction

**Hmmix-nf** is a nextflow implementation of [hmmix](https://github.com/LauritsSkov/Introgression-detection/tree/master), a method for identifying archaic introgression in modern human genomes. 

## Pipeline Summary

- **Input:** VCF/BCF files, callability mask BED file, ancestral allele FASTA, reference genome FASTA, sample JSON.

**1.** Find derived variants in outgroup (can be skipped by using pre-made file available [here](https://zenodo.org/records/11212339)).

**2.** Estimate local mutation rate (can be skipped by using pre-made file available [here](https://zenodo.org/records/11212339)).

**3.** Find variants in ingroup.

**4.** Train the HMM.

**5.** Decoding.
- **Output:** .txt file for each individual with a list of segments annotated to one of two states: Human or Archaic

