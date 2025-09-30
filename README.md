# Transcriptomics Casus: Genexpressieanalyse bij Reumatoïde Artritis

<p align="center">
  <img src="assets/0419_RA-Symptoms.jpg" alt="Reuma Symptomen" width="600"/>
</p>

## Inhoud/structuur

- `Bam bestanden/` - reads dat gemapt zijn tegen referentiegenoom  
- `bronnen/` - gebruikte bronnen 
- `RScripts/` – Rscript waar alle data wordt verwerkt 
- `resultaten/` - grafieken en tabellen
- `Ruwe Data/` – Hierin zit de data van de Patiënten met en zonder Reumatoïde Atritis (RA).
- `assets/` - overige documenten voor de opmaak van deze pagina
- `README.md` - Is het document waar het verslag (de tekst) in wordt gemaakt. 






---

## Inleiding

Reumatoïde artritis (RA) is een chronische, systematische auto-immuunziekte waarbij het immuunsysteem het eigen lichaam aanvalt. De exacte oorzaak is nog onbekend, maar vermoedelijk speelt een combinatie van genetische aanleg, omgevingsfactoren en een ontspoord immuunsysteem een rol [Gabriel, 2001](Bronnen/gabriel2001.pdf). Een belangrijk kenmerk van RA is synovitis: een ontsteking van het gewrichtsslijmvlies, wat leidt tot pijn, zwelling en uiteindelijk gewrichtsschade [Radu & Bangua, 2021](Bronnen/cells-10-02857-v2.pdf). Vroege opsporing en behandeling kunnen schade beperken, maar genezing is (nog) niet mogelijk.

In dit project is RNA-sequencing toegepast op synoviumbiopten van zowel gezonde personen als patiënten met RA (gevestigde diagnose >12 maanden), met als doel het identificeren van verschillen in genexpressie en betrokken biologische processen. De analyse is uitgevoerd in [Rstudio](RScripts/ProjectRA.R) en omvat onder meer differentiële expressieanalyse en functionele verrijkingsanalyse op basis van Gene Ontology (GO) en KEGG-pathways.

De gebruikte brondata en artikelen zijn te vinden in de folder [bronnen](Bronnen). 



## Methoden

Voor deze studie is RNA-sequencing data geanalyseerd van synoviumbiopten van gezonde controles en patiënten met reumatoïde artritis (RA).

<p align="center">
  <img src="assets/BiorenderFlowchart2.png" alt="flowchart" width="600"/>
</p>

_Figuur 1: flowschema van dataverwerking in R._


#### Verkregen data

De RNA-sequencing data zijn verkregen uit synoviumbiopten van 4 patiënten met RA en 4 gezonde controles. Bij de RA-patiënten is de diagnose bevestigd door de aanwezigheid van anti-CCP autoantistoffen. Een overzicht van de data is te vinden in  ruwe data link(data link ruwe data).

_Tabel 1. overzicht van de gebruikte samplemetadate in deze studie._

| Sample ID     | Age | Sex    | Condition                          |
|---------------|-----|--------|------------------------------------|
| SRR4785819    | 31  | Female | Normal                             |
| SRR4785820    | 15  | Female | Normal                             |
| SRR4785828    | 31  | Female | Normal                             |
| SRR4785831    | 42  | Female | Normal                             |
| *Average Age* | *29.8* | – | *Normal*                       |
| SRR4785979    | 54  | Female | Rheumatoid arthritis (established) |
| SRR4785980    | 66  | Female | Rheumatoid arthritis (established) |
| SRR4785986    | 60  | Female | Rheumatoid arthritis (established) |
| SRR4785988    | 59  | Female | Rheumatoid arthritis (established) |
| *Average Age* | *59.8* | – | *Rheumatoid arthritis (established)* |


#### Mappen van data en countmatrix 
De ruwe sequencingbestanden (FASTQ-formaat) zijn opgeslagen in de map [Ruwe Data](Ruwe%20Data). Voor de uitlijning is het humane referentiegenoom GRCh38 (GCF_000001405.40, versie GRCh38.p14) van NCBI gebruikt. Vanwege de grootte van het genoom is het FASTA-bestand niet opgenomen in de repository, maar kan via NCBI worden gedownload via accessionnummer [GCF_000001405.40](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/).

De referentie-index is opgebouwd met het R-pakket Rsubread (versie 2.20.0) met voldoende geheugen (4 GB) om de indexbestanden te splitsen. De sequencing reads zijn uitgelijnd met de functie align(), waarna de resulterende BAM-bestanden gesorteerd en geïndexeerd zijn met Rsamtools. De genexpressie is geteld met featureCounts(), waarbij gebruik is gemaakt van de GTF-annotatie Homo_sapiens.GRCh38.114.chr.gtf.gz van ENSEMBL.

Differentiële expressieanalyse is uitgevoerd met DESeq2, waarbij het model ~ treatment (controle vs. RA) is toegepast. Resultaten zijn gevisualiseerd met EnhancedVolcano en geanalyseerd met goseq en pathview om biologische functies en pathways te identificeren.

#### Statistische analyse
De gegenereerde countmatrix (data/count_matrix.txt) en een behandel-tabel met controle- en RA-status zijn ingeladen in DESeq2 (versie 1.46.0) om differentiële genexpressie te berekenen. Hierbij zijn log₂ fold changes, p-waarden en meervoudige testcorrecties (Benjamini-Hochberg) bepaald. Resultaten zijn gevisualiseerd in een volcano plot (log₂ fold change vs. significatie). Daarnaast is functionele verrijkingsanalyse uitgevoerd met goseq voor Gene Ontology en met pathview voor KEGG pathway-analyse.

## Resultaten

### Differentiële genexpressieanalyse 
om genen te identificeren die verschillend tot expressie kwamen tussen RA-patiënten en gezonde controles, werd een RNA-seq analyse uitgevoerd met DESeq2. In totaal werden 102 genen met verhoogde expressie (log₂FC > 1, padj < 0.05) en 88 genen met verlaagde expressie (log₂FC < -1, padj < 0.05) geïdentificeerd. de resultaten worden weergegeven in een [volcano](Resultaten/Plots/VolcanoplotRA.png) , waarin de log2 fold change tegenover de -log1- p-waarde van elk gen staat. De volledige lijst van significante genen is beschikbaar in [ResultatenRA](Resultaten/ResultatenRA.csv).

De plot toonde aan dat er meerdere genen significant op- of neerwaarts gereguleerd zijn. echter ontbraken labels voor de belangrijkste genen, waardoor het lastig was om direct te zien welke genen betrokken waren bij RA. voor toekomstige analyse zou het toevoegen van genlabels helpen om de interpretatie te verbeteren. 
Door technische problemen in R kon dit op dit moment niet worden gerealiseerd. 
 
### GO-enrichmentanalyse
Om de biologische betekenis van de differentieel tot expressie gebrachte genen te achterhalen, werd een [GO-enrichmentanalyse](Resultaten/Plots/Rplot03.png) uitgevoerd. De resultaten werden weergegeven in een barplot waarin termen zoals: _Immune system process, Immune respons, Protein binding, Intracellular orgaanelle lumen_ 
Duidelijk naar voren komen. deze termen waren statistisch verrijkt en wezen op verstoring van immuunprocessen, wat goed paste bij het ziektebeeld van RA. 
De visualisatie toonde zowel het aantal betrokken genen als p-waarden, en gaf een helder overzicht van de belangrijkste biologische processen.

### KEGG-pathwayanalyse
voor verdere interpretatie werd een [KEGG-pathwayanalyse](Resultaten/hsa03260.png) uitgevoerd. het bijgevoegde figuur toonde een HIV-gerelateerde pathway, wat niet direct relevant leek voor RA. Toch bevatte deze pathwat immuunreceptoren zoals _CD4_, die ook een rol spelen in auto-immuunziekten zoals RA. 

Hoewel deze pathway enkele relevante elementen bevat, zou het beter zijn om pathways te tonen die direct met RA te maken hebben, zoals: 
_Cytokine-cytokine receptor interaction, Toll-like receptor signaling pathway_ , deze zijn nauw betrokken bij ontstekingsreacties en geven een duidelijker beeld van de moleculaire processen in RA.


## Conclusie
Deze RNA-sequencinganalyse van synoviumbiopten van RA-patiënten en gezonde controles laat zien dat er duidelijke verschillen zijn in genexpressie tussen beide groepen. De differentiële expressieanalyse identificeerde meerdere genen met significante op- of neerregulatie. hoewel verdere interpretatie beperkt werd door het ontbreken van genlabels in de visualisatie.

De GO-enrichmentanalyse toonde een duidelijke verrijking van immuun-gerelateerde processen, zoals _immune respone_ en _Protein binding_, wat goed aansluit bij het ontstekingskarakter van RA. Dit bevestigd dat immuunactivatie een centrale rol speelt in pathogenese van de ziekte.

Hoewel de KEGG-pathwayanalyse oorspronkelijk een HIV-gerelateerde pathway toonde, bevatte deze wel immuunreceptoren die ook relevant zijn voor RA. voor een meer directe interpretatie zouden RA-specifieke pathwats zoals _Cytokine-cytokine receptor interaction_ beter geschikt zijn. 

Afsluitend bieden de resultaten waardevolle inzichten in de moleculaire processen die betrokken zijn bij RA, met name op het gebied van immuunregulatie. verdere verfijning van de visualisaties en pathwayselectie kan bijdragen aan een nog scherpere interpretatie en mogelijke aanknopingspunten voor therapie of diagnostiek.




