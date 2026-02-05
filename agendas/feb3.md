# Team 7 Meeting #4, February 3rd, 2026
Research Question:
> Do fish species and their associated microbial diversity in hindgut samples influence swimming performance and preferred depth in water? 

##  Fish Dataset and Group Assignment 2 - Minutes by Izaak Yip

### Develop 'Aims'
- Focused questions to answer specific parts of the main research question

#### Aims
1. Alpha & beta diversity correlation between microbiome & performance
2. Core microbiome analysis, Different locations have different microbiota
3. Functional analysis - picrust2
-  M1: see pathways are upregulated, check pathways to see if they produce more energy - MetaCyc or Kegg
-  M2: going in blind or literature search first
4. DEseq (differential abundance) - inform the relations to picrust (If one taxa is upregulated or not)
- Comparison of singular groups between different water levels - subgroup data - trophic levels (predators vs prey fish), breakdown

##  Comments from Meeting
For (1) Alpha & beta diversity metrics
Can break them down further into...
- trophic levels
- substrata collection (although this may become less relevant, since seafloor is a tangent to our research question)
For (2) Compositional metrics on core microbiome
- other students found different geographical locations had different specific core microbiomes
For (3) Functional analysis
- Consider whether fast fish have certain microbiota
- which pathways are regulated
- for literature search, check which pathways to look for (KEGG, etc.)
- for KEGG, a good starting point would be to check pathways associated with higher energy outputs
- e.g. for humans there's SCFA production, complex carb breakdown
##### Subgrouping our data
Imogen suggested another possible addition to our work once completing initial analyses
If one group (e.g. fast high fish) stands out in terms of significantly higher diversity, we can run specific analyses on that group.
Within this group... how does a certain variable (decide on any variable) affect them?

# Action steps
1/ check correlation between speed and depth (e.g. fast high fish, fast low fish, etc.)
2/ complete module on functional analysis


_note last edited Feb 5 by VX_


### Running Analyses
- Which analyses would be the most appropriate to run on this dataset, given the read quality and sample size (from Group Assignment 1)

Two Options
1. Find and know the pathways we need to analyze.
- Metasiq or keg and see what patterns exist for metabolites
- Conduct a literature search
2. Blind do analysis and see patterns emerge

- Make a comparison table of different variables

### Background for Fish
- Given that this is a paper meant for other scholars, how much background information is expected for this manuscript?
- Discuss fish taxa?
- Importance of fish in our ecosystem?

**Comments**\
The research papers cited in creating the research questions can be found here.
- https://pubmed.ncbi.nlm.nih.gov/39466622/
- https://pubmed.ncbi.nlm.nih.gov/40564368/
### Action Items
|Person|Task       |Done?|
|:---: |:---:      |:---:|
|IY    |Code Fix   |[ ]  |
|IY    |`ggpairs`  |[ ]  |
|SD    |Graphs     |[ ]  |
|AL    |Diversity  |[ ]  |
|VH    |Lit. Review|[ ]  |
|VX    |Lit. Review|[ ]  |

