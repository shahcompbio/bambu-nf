# shahcompbio/bambu-nf: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 1.1.0 - [2025-07-27]

Updates to perform QC post-filtering; make filtering the default pre-assembly

### `Added`

- Perform QC with nanostats post-filtering
- filter reads pre-assembly as default behavior
- quant only mode
- merge only mode

## 1.0.0 - [2025-05-27]

Initial release of shahcompbio/bambu-nf, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Transcript assembly & quant with bambu
- Flags to pre-filter reads before assembly on read length & mapq (useful for samples with inadequate QC)
- Flags to pre-filter reads on accessory chromosomes (can sometimes cause issues for Bambu)
- Flag to adjust NDR in bambu
