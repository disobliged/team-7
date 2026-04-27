# Exploring Hindgut Microbiome Composition and Predicted Energy Metabolism Pathways Across Fish with Contrasting Swimming Performance
![RStudio](https://img.shields.io/badge/RStudio-4285F4?style=for-the-badge&logo=rstudio&logoColor=white)
![R](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![sponsors](https://img.shields.io/opencollective/sponsors/0?color=t&label=Generous%20Sponsors)

***

Analyses conducted in R to research how the hindgut microbiome affects fish swim speed. This project was made for MICB 305 and has been published on \[UJEMI].

## Research Question
> Does hindgut microbial diversity, core microbiota, differentially abundant taxa and predicted energy metabolic pathways correlate with swimming performance across fish species?

## Contributors!
<a href="https://github.com/disobliged/Fish-Swim-Performance/graphs/contributors">
  <img src="https://contrib.rocks/image?repo=disobliged/Fish-Swim-Performance" />
</a>

*Made with [contrib.rocks](https://contrib.rocks).*

### Mentors
We received generous support from [Avril Metcalfe-Roach](https://github.com/armetcal) and [Imogen Porter](https://github.com/iporter-16) in the creation of our project. We thank them for their support during the term.

## Agenda Formatting
Both agendas and meeting minutes are consolidated into the same file for efficiency and tidiness.

## GitHub Directory

```
.
├── agendas
├── datasets
│   ├── picrust_output
│   │   ├── EC_metagenome_out
│   │   ├── intermediate
│   │   │   ├── pathways
│   │   │   │   └── minpath_running
│   │   │   ├── place_seqs_arc
│   │   │   │   └── epa_out
│   │   │   └── place_seqs_bac
│   │   │       └── epa_out
│   │   ├── KO_metagenome_out
│   │   └── pathways_out
│   └── picrust_strat_output
│       ├── EC_metagenome_out
│       ├── intermediate
│       │   ├── pathways
│       │   │   └── minpath_running
│       │   ├── place_seqs_arc
│       │   │   └── epa_out
│       │   └── place_seqs_bac
│       │       └── epa_out
│       ├── KO_metagenome_out
│       └── pathways_out
├── plots
│   └── archive
└── scripts
    └── archivedScripts
```

- ```agendas``` is where our meeting agendas/minutes live.
- ```datasets``` is where we stored the datasets we used and RDS files we made.
- ```picrust_output``` and ```picrust_strat_output``` is where the PICRUSt2 unstratified and stratified outputs are saved respectively.
- ```plots``` is where our generated plots are saved.
- ```scripts``` is where we made our scripts.
- The ```archive``` and ```archivedScripts``` folders are where we saved unused plots and scripts respectively

## Open Source
The project's dataset is open access and can be found [here](https://doi.org/10.1038/s41467-022-34557-2).
In an effort to keep research accessible, the project's code is licensed under the GPL v3 license. This is done so the code will always be free and accessible for future reference.
