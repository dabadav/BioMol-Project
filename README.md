# Molecular Biology

### Computational analysis Project - Group 11

Dante Aviñó, Miguel Borge, Pol Segura

------

The main objective of the project is to **annotate the genes** included in the contig and to **functionally characterize** the predicted proteins.

Download assigned contig from a nematode species.

```bash
# Using npm package manager we're going to install github-files-fetcher
# so we can download particular files from the github repo.
npm install -g github-files-fetcher

# Download contig sequence from repo
fetcher --url="https://github.com/dantekali/BioMol-Project/blob/main/Group11_contig_194888_195063.fa"  --out="~/Desktop/Project"
```



------

#### Genes annotation included in the contig

We will first perform a **blastx** of the contig, against the *nr database* to get some insight of the possible *gene* it might code, and the *species* it belongs.
  - https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastx&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome

> We download the aa sequence of the first protein match (754 aa) into a fasta format file called "[CAEBREN_29266.fa](https://github.com/dantekali/BioMol-Project/blob/main/CAEBREN_29266.fa)"

| \>EGT37340.1 hypothetical protein CAEBREN_29266 [Caenorhabditis brenneri], 754 aa |
| :----------------------------------------------------------- |
| MGRYVDMLSSTSPLEFTKVIVASLDYSNEGMTRVILRKALTSASESSRKWTTRYLAVLASSDLPMFSDWGIQLMLRQLADESSKVVRHTIRILSRWLPEHPSRNLRKCEWSVFGEAGDLLKAHVYALEFECASDEDEVRDVIRFWMTDFNKKYLQIIDEEMKEMMFHVKRSIDGSFSRSSSDRPDTSLGVHAPLHLFAALGGHETGKRILLEENVCEELLSVIRIGKCFEELKSSLLALASIGSTDRGFEILPLDAVPTVLKIAEEHTVLTVRGIAFWALCTFSQCIEGAKRLAAFGWESNRFRYAMDIARGKISEDEGMISTPVAGTTAGSVSSTWRPARKITMQHHRHSSLFDSQINVKQSRAKSESAVSRRGNSKGRRRSQSEGDIQEKSPKRESRIDSFFSQRLWNSEKYLYKSSGTSDSSSITYHKRTVTNSSSGYHIQEEITVTVSPPGHLFPDESVAKSAATSRLSTDRRRANTTNSLFEEEEAPKTRSSTVARCIREGLKITSEELEAEGVVADTIMEPHFSCRLREKYHLMPFRVRACLHINRHVGDPIRYVFMTREEERHFADYRRQVLHDPWLFNELRKEDNAVKKTINVVPLQTVALPTEIEIMCGNIFPAKPKSDPIFSFHENDDSAVEDRGARTGHARSGIHIQPHSAYRCFHCSSNEDSVRGYPHPDAPMLRKEVLGQVDMLEIKEYPAKRLIGLRQHNPWLFQWPCMYADVLELLDEYRFKPHSRAFLQHIFYDALQI |

| Select for downloading or viewing reports |                         Description                          |       Scientific Name       | Max Score | Total Score | Query Cover | E value | Per. Ident | Acc. Length | Accession  |
| :---------------------------------------: | :----------------------------------------------------------: | :-------------------------: | :-------: | :---------: | :---------: | :-----: | :--------: | :---------: | :--------: |
|        Select seq gb\|EGT37340.1\|        | hypothetical protein **CAEBREN_29266** [Caenorhabditis brenneri] | **Caenorhabditis brenneri** |   1335    |    1413     |     24%     |   0.0   |   86.87%   |     754     | EGT37340.1 |

We will obtain gene predictions using the two types of methods and we will compare the results with Ensembl annotations.

1. **ab-initio tools** to obtain a first prediction (select closest species in model)
    <!--Ab-initio methods: they use several elements in the genomic sequence (suchas donor and acceptor splice sites, branch site, initiation and termination codons)and codon usage to obtain a model based on a training set.-->

  - <u>GeneID  prediction</u>: https://genome.crg.cat/software/geneid/geneid.htmlWe 

    We use the "contig.fa" file to make a gene prediction using GeneID, this will result in the generation of a gff file called "[GeneID.gff](https://github.com/dantekali/BioMol-Project/blob/main/GeneID.gff)" containing in this case 4 predicted proteins of different aa length.
    
| >GeneID\|Protein 1\| 2 exons, 56_AA [Forward]                |
| :----------------------------------------------------------- |
| GFLDYSYGNMFGANDESFYSRLLMSLPSIIFSFLLNEKGKADEMIEWRTSSSSGG*     |
| **>GeneID\|Protein 2\|  4 exons, 329_AA **[Reverse]          |
| MNKKTKYVLGGVVFSSSIAIVAAIFACGFMVFDISEFQNNIRSDLKEFQFYSSDSWNVMLKGKTIGFRVRRQYPSAPVGGGGQIEDGTCQCAEQSTGCPPGPPGPSGTPGHPGDSGAPGNPGQPGSAGIVEMHESMKNGCISCPQGPAGSPGPDGPPGPPGPSGNPGRESPAGPAGQPGPPGGLGPPGQNGNPGAPGNMGAPGKPGMKHTNPPGQPGPTGPMGPPGPPGNDAQFANGPPGPPGPMGPPGKPGSAGKDGQDGNPGSDGHPGSDGQYCPCPSRTPNLGVNGFSQDEENAVDTSFGFTKRNLVKMKRMMAKLHKKFSSAIA* |
| **>GeneID\|Protein 3\| 7 exons, 784_AA **[Reverse]           |
| MPRLVRDDDDSRNNRLHLSLCVLWSTAFLLLRYVDMLSSTSPLEFTKVIVASLDYSNEGMTRVILRKALTSASESSRKWTTRYLAVLASSDLPMFSDWGIQLMLRQLADESSKVVRHTIRILSRWLPEHPSRNLRKCEWSVFGEAGDLLKAHVYALEFECASDEDEVRDVIRFWMTDFNKKYLQIIDEEMKEMMFHVKRSIDGSFSRSSSDRPDTSLGVHAPLHLFAALGGHETGKRILLEENVCEELLSVIRIGKCFEELKSSLLALASIGSTDRGFEILPLDAVPTVLKIAEEHTVLTVRGIAFWALCTFSQCIEGAKRLAAFGWESNRFRYAMDIARGKISEDEGMISTPVAGTTAGSVSSTWRPARKITMQHHRHSSLFDSQINVKQSRAKSESAVSRRGNSKGRRRSQSEGDIQEKSPKRESRIDSFFSQRLWNSEKYLYKSSGTSDSSSITYHKRTVTNSSSGYHIQEEITVTVSPPGHLFPDESVAKSAATSRLSTDRRRANTTNSLFEEEEAPKTRSSTVARCIREGLKITSEELEAEGVVADTIMEPHFSCRLREKYHLMPFRVRACLHINRHVGDPIRYVFMTREEERHFADYRRQVLHDPWLFNELRKEDNAVKKTINVVPLQTVALPTEIEIMCGNIFPAKPKSDPIFSFHENDDSAVEDRGARTGHARSGIHIQPHSAYRCFHCSSNEDSVRGYPHPDAPMLRKEVLGQVDMLEIKEYPAKRLIGLRQHNPWLFQWPCMYADVLELLDEYRFKPHSRAFLQHIFYDALQI* |
| **>GeneID\|Protein 4\| 3 exons, 156_AA** [Forward]           |
| MSVSTIRAISILLLLASYTLADLPSCARAKCVHCAVDFIDRMCPTACAGCKTTHQSIQHCTYMLQAVNIQRPPQPPPAFNTQSQFTQSVNTRQISNEGGPQVQRPQPQPVHPQPQQPQQQ21HQQQFQQQQQFQPQVQTQQQPQQPLLNIPLQPHAQP |

  - <u>FGENESH</u>: http://www.softberry.com/berry.phtml?topic=fgenesh&group=programs&subgroup=gfindYou

    We use the "contig.fa" file to make a gene prediction using FGENESH, this will result in the generation of a txt file called "[FGENESH.txt](https://github.com/dantekali/BioMol-Project/blob/main/FGENESH.txt)", containing in this case 3 predicted proteins of different aa length.

| \>FGENESH:   1   3 exon (s)   1526  -   2189   161 aa, chain + |
| :----------------------------------------------------------- |
| AAENSFSPSNFVCIFKRSLFEKALNCLAKTEPITFLKCDHECHEEIVRSDRLKQVTPNIQNQIFISSELATYETELDKLCRFQSCYMNCMAPVVKEMCGEDESKHAVEIVESYVQWHADDISDWHSITGNDETLPKACQTLVKTHSKADDPILQIIGNAAL |
| **\>FGENESH:   2   7 exon (s)   2396  -   6562   722 aa, chain -** |
| MGRYVDMLSSTSPLEFTKVIVASLDYSNEGMTRVILRKALTSASESSRKWTTRYLAVLASSDLPMFSDWGIQLMLRQLADESSKVVRHTIRILSRWLPEHPSRNLRKCEWSVFGEAGDLLKAHVYALEFECASDEDEVRDVIRFWMTDFNKKYLQIIDEEMKEMMFHVKRSIDGSFSRSSSDRPDTSLGVHAPLHLFAALGGHETGKRILLEENVCEELLSVIRIGKCFEELKSSLLALASIGSTDRGFEILPLDAVPTVLKIAEEHTVLTVRGIAFWALCTFSQCIEGAKRLAAFGWESNRFRYAMDIARGKISEDEGMISTPVAGTTAGSVSSTWRPARKITMQHHRHSSLFDERHFADYRRQVLHDPWLFNELRKEDNAVKKTINVVPLQTVALPTEIEIIPRLPTSRCSNASKRGSWSSGYVRNQRVSSKEIDWVRDYQFYSSDSWNVMLKGKTIGFRVRRQYPSAPVGGGGQIEDGTCQCAEQSTGCPPGPPGPSGTPGHPGDSGAPGNPGQPGSAGIVEMHESMKNGCISCPQGPAGSPGPDGPPGPPGPSGNPGRESPAGPAGQPGPPGGLGPPGQNGNPGAPGNMGAPGKPGMKHTNPPGQPGPTGPMGPPGPPGNDAQFANGPPGPPGPMGPPGKPGSAGKDGQDGNPGSDGHPGSDGQYCPCPSRTPNLGVNGFSQDEENAVDTSFGFTKRNLVKMKRMMAKLHKKFSSAIA |
| **>FGENESH:   3   3 exon (s)   9735  -  10306   108 aa, chain +** |
| MCPTACAGCKTTHQSIQAVNIQRPPQPPPAFNTQSQFTQSVNTRQISNEGGPQVQRPQPQPVHPQPQQPQQQHQQQFQQQQQFQPQVQTQQQPQQPLLNIPLQPHAQP |

  - <u>GENESCAN</u>: http://argonaute.mit.edu/GENSCAN.html

    We use the "contig.fa" file to make a gene prediction using GENESCAN, this will result in the generation of a txt file called "[GENESCAN.txt](https://github.com/dantekali/BioMol-Project/blob/main/GENESCAN.txt)".

| >GENSCAN_predicted_peptide_1\|132_aa                         |
| :----------------------------------------------------------- |
| TEPITFLKCDHECHEEIVRSDRLKQVTPNIQNQIFISSELATYETELDKLCRFQSCYMNCMAPVVKEMCGEDESKHAVEIVESYVQWHADDISDWHSITGNDETLPKACQTLVKTHSKADDPILQIIGNAAL |
| **>GENSCAN_predicted_peptide_2\|962_aa**                     |
| XCQKQQNADGTYRRNRHRSLKFQRSMKKKSVCGRNSQRVWIKRERRYVDMLSSTSPLEFTKVIVASLDYSNEGMTRVILRKALTSASESSRKWTTRYLAVLASSDLPMFSDWGIQLMLRQLADESSKVVRHTIRILSRWLPEHPSRNLRKCEWSVFGEAGDLLKAHVYALEFECASDEDEVRDVIRFWMTDFNKKYLQIIDEEMKEMMFHVKRSIDGSFSRSSSDRPDTSLGVHAPLHLFAALGGHETGKRILLEENVCEELLSVIRIGKCFEELKSSLLALASIGSTDRGFEILPLDAVPTVLKIAEEHTVLTVRGIAFWALCTFSQCIEGAKRLAAFGWESNRFRYAMDIARGKISEDEGMISTPVAGTTAGSVSSTWRPARKITMQHHRHSSLFDSQINVKQSRAKSESAVSRRGNSKGRRRSQSEGDIQEKSPKRESRIDSFFSQRLWNSEKYLYKSSGTSDSSSITYHKRTVTNSSSGYHIQEEITVTVSPPGHLFPDESVAKSAATSRLSTDRRRANTTNSLFEEEEAPKTRSSTVARCIREGLKITSEELEAEGVVADTIMEPHFSCRLREKYHLMPFRVRACLHINRHVGDPIRYVFMTREEERHFADYRRQVLHDPWLFNELRKEDNAVKKTINVVPLQTVALPTEIEIMCGNIFPAKPKSDPIFSFHENDDSAVEDRGARTGHARSGIHIQPHSAYRCFHCSSNEDSGKTIGFRVRRQYPSAPVGGGGQIEDGTCQCAEQSTGCPPGPPGPSGTPGHPGDSGAPGNPGQPGSAGIVEMHESMKNGCISCPQGPAGSPGPDGPPGPPGPSGNPGRESPAGPAGQPGPPGGLGPPGQNGNPGAPGNMGAPGKPGMKHTNPPGQPGPTGPMGPPGPPGNDAQFANGPPGPPGPMGEFLRCDRVYVIIIMFSGPPGKPGSAGKDGQDGNPGSDGHPGSDGQYCPCPSRTPNLGVNG |



<u>Sequences Comparison</u>: compare them with the annotated protein in found in blastx nr database: [CAEBREN_29266.fa](https://github.com/dantekali/BioMol-Project/blob/main/CAEBREN_29266.fa)

To obtain the multiple alignment we will use the simple MSA option from T-coffee software (http://tcoffee.crg.cat/).

We first observe the results of the different software looking for proteins coded in the same region of contig with similar polypeptide chain lengths.

">CAEBREN_29266 [Caenorhabditis brenneri], 754 aa"
">GeneID|Protein 3| 7 exons, 784_AA [Reverse]"
**">FGENESH:   2   7 exon (s)   2396  -   6562   722 aa, chain -"**
">GENSCAN_predicted_peptide_2\|962_aa"

As we can see the putative proteins: GeneID_Protein3, FGeneSH_Protein2 and the blastx match CAEBREN_29266 protein have similar lengths and same exon number. 

In this case GENSCAN didn't give us useful prediction. So we will use the T-Coffee MSA tool to compare the selected predictions from GeneID and FGENESH against the CAEBREN blasted protein.

- We first compare the GeneID protein with the CAEBREN_29266 one, we put them in the same text file: [tcoffee_geneid.txt](https://github.com/dantekali/BioMol-Project/blob/main/tcoffee_geneid.txt)
- Then FGENESH protein with the CAEBREN_29266. We create the t-coffee input file by putting them in the same text file: [tcoffee_fgenesh.txt](https://github.com/dantekali/BioMol-Project/blob/main/tcoffee_fgenesh.txt)

By analyzing the [FGENESH](https://github.com/dantekali/BioMol-Project/blob/main/result_fgenesh.pdf) and [GeneID](https://github.com/dantekali/BioMol-Project/blob/main/result_geneid.pdf) results we conclude that the protein_2 predicted by FGENESH is almost identical to the blast result, in the other hand GeneID protein prediction presents more variations.




2. **homology-based tools** with the annotations of  a closely related species (web-server/local)
    <!--Gene predictions are based on alignments from known proteins (usually) from other genomes.-->

  - Run *blastx* of the unspliced sequence against the *nr database*. Since this is a huge database, this time we will run blast on the ncbi server https://blast.ncbi.nlm.nih.gov/Blast.cgi
  - https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome
  > We transform the translated nucleotide sequence "contig.fa" to protein with blastx, in this step we generate a "tioblast.fa" containing the putative proteins.
```bash
bedtools getfasta -fi Code.fa -bed GeneID.gff>geneID.fa

exonerate     -m     p2g     --showtargetgff     -q     tioblast.fa     -t Code.fa

exonerate     -m     p2g     --showtargetgff     -q     tioblast.fa     -t Code.fa -S F

exonerate     -m     p2g     --showtargetgff     -q     tioblast.fa     -t Code.fa -S F| egrep -w exon > GeneID.gff

bedtools getfasta -fi Code.fa -bed GeneID.gff>exonerate.fa

sed     -e    '2,$s/>.*//'    exonerate.fa     | grep     -v    '^$'     >exonerate_singleLine.fa
```
  > We transform "tioblast.fa" file containing the putative proteins with exonerate finally outputing a exonerate_singleLine.fa



------

#### Functionally characterization of the predicted proteins

Functional annotation is an essential step in omics data analysis. It is defined as the process of attaching biological information to gene and protein sequences (such as those predicted using Ab-initio and Homology-based methods) and describing functional groups based on similarities.


- Search the Gene Ontology database (http://geneontology.org)
- Identify conserved protein domains in the peptides encoded by genes (you can find the FASTA file with the sequences -query.fas- here. This python script allows the programmatic access to InterPro Web Services using an API (Application Programming Interface). Use the script  iproscan5.py
- Use   the   script  find_enrichment.py (from the  goatools python library;https://github.com/tanghaibao/goatools) to find if the genes under study (e.g. genes expressed in a specific tissue)  are enriched in a particular (or various) GO terms. For a GO enrichment analysis you need to know 
- i) the names or database identifiers of the genes you want to test(in this example the InterPro domains found for all genes specifically expressed in the tissue analyzed in our study -study.names), 
- ii) the names or database identifiers of the genes in the population(e.g. the InterPro domains of all genes expressed in this organism - i.e., in all tissues - population.names), 
- iii) an association file that maps a gene or database identifier name to a GO category (InterPro2GO.map) and 
- iv) a file with a basic version of the GO database (go-basic.obo)



#### Discuss the performance of the different methods. Can you guess to which species does your contig correspond?
