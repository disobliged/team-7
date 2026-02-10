# Team 7 Meeting #5, February 10th, 2026
Research Question:
> Do fish species and their associated microbial diversity in hindgut samples influence swimming performance and preferred depth in water?


### Action Items Recap
Person|Task       |Done?|
|:---: |:---:      |:---:|
|IY    |Code Fix   |[✅]  |
|IY    |`ggpairs`  |[✅]  |
|SD    |Graphs     |[✅]  |
|AL    |Diversity  |[✅]  |
|VH    |Lit. Review|[✅]  |
|VX    |Lit. Review|[✅]  |



## Updated Plot:

<img width="526" height="450" alt="Screenshot 2026-02-09 at 10 36 54 PM" src="https://github.com/user-attachments/assets/21520675-653e-48c0-b3dd-b03b20656233" />

### Notes: 
- Make sure to pick a number that maintains a certain amount of samples, that doesn't lose a lot of features (maximizing).


## Depth by Swim Performance (boxplot):

<img width="526" height="450" alt="Screenshot 2026-02-09 at 10 40 06 PM" src="https://github.com/user-attachments/assets/4bf33221-51b4-4603-8f3c-1772cb9a304b" />

### Notes: 
- Seems like there is a trend and have their own biases.
- Because a lot of them have a trend, we will have a spread and group them into shallow and deep.
- A lot of fish are generalists or caught a bunch of generalists.
- Discussion: normalizing by species - bias towards more common fish.

## Performance Counts: Shallow vs Deep (30m+) (bar plot):

<img width="526" height="450" alt="Screenshot 2026-02-09 at 10 41 40 PM" src="https://github.com/user-attachments/assets/81f7bd79-9ea3-4a03-8883-70ef2b093d42" />

## Literature Search 

MetaCyc vs KEGG
- Most studies seem to use KEGG rather than MetaCyc
- In the context of our research question, MetaCyc might be more suitable because it provides detailed information on individual enzyme reactions and steps of specific metabolic pathways (rather than the more general trends that KEGG would provide)
- The pathways in KEGG are standardized and we could still determine which ones affect swimming performance and preferred depth in water but we would not have a layered understanding of the specific steps/enzymes in the metabolic processes which drive these outcomes

### Notes: 
- KEGG is outside of metabolism
- Always start with KEGG, then work pathways are at a broader scale. Looking for specific high-energy metabolic pathways, that are relevant.
- Step back and look at MetaCyc, the pathway may not be all significantly upregulated (harder to find pathways vs. steps)

## Functional Analysis Approach
1. Find and know the pathways we need to analyze
   - Pros: Perhaps easier to interpret and understand patterns because they already established in the literature, may strenghten existing findings
   - Cons: A more bias approach which limits the exploration of new patterns 
  
2. Going in blind (without literature guidance)
   - Pros: Might find new and novel patterns, allowing trends to emerge naturally without inherent bias from previous literature, minimizes bias
   - Cons: Might find biologically irrelevant patterns, trends might be harder to interpret without context
  
   - -> Decide which approach to use
   - PICRUSt good for 16S based functional prediction

Zoom link: https://ubc.zoom.us/j/9684360467?pwd=bU5ORnh0Z0I4QkZFVkVkVWRDbDhnUT09

### Notes: 
- Physicality and movement, we have some idea of what pathways we would analyze. Literature dive and look for specific pathways.
- Blind approach (top 50 pathways and see which ones are of interest), or just the interesting ones.
- Go in informed, then blind
- Clear idea of what we hope to see
- PICRUSt would be good
- Look at modules beforehand, look ahead and know what's happening

## Group 2 Assignment - Analyses

Fill in analysis table:

For each technique, decide:
- a) Which ones you currently plan to use (highlight, change colour, etc)
- b) Whether each chosen technique typically uses statistics*
- c) Which aim each chosen technique will contribute to
- d) How each chosen technique will help you address your research question biologically
