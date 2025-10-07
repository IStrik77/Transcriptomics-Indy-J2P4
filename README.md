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

Reumatoïde artritis (RA) is een chronische, systemische auto-immuunziekte waarbij het immuunsysteem het eigen lichaam aanvalt. De exacte oorzaak is nog onbekend, maar vermoedelijk speelt een combinatie van genetische aanleg, omgevingsfactoren en een ontspoord immuunsysteem een rol [Gabriel, 2001](Bronnen/gabriel2001.pdf). Een belangrijk kenmerk van RA is synovitis: een ontsteking van het gewrichtsslijmvlies, wat leidt tot pijn, zwelling en uiteindelijk gewrichtsschade [Radu & Bangua, 2021](Bronnen/cells-10-02857-v2.pdf). Vroege opsporing en behandeling kunnen schade beperken, maar genezing is (nog) niet mogelijk.

In dit project is RNA-sequencing toegepast op synoviumbiopten van zowel gezonde personen als patiënten met RA (vastgestelde diagnose >12 maanden), met als doel het identificeren van verschillen in genexpressie en betrokken biologische processen. De analyse is uitgevoerd in [Rstudio](RScripts/ProjectRA.R) en omvat onder meer differentiële expressieanalyse en functionele verrijkingsanalyse op basis van Gene Ontology (GO) en KEGG-pathways.

De gebruikte brondata en artikelen zijn te vinden in de folder [bronnen](Bronnen). 



## Methoden

Voor deze studie is RNA-sequencing data geanalyseerd van synoviumbiopten van gezonde controles en patiënten met reumatoïde artritis (RA).

<p align="center">
  <img src="assets/BiorenderFlowchart2.png" alt="flowchart" width="600"/>
</p>

_Figuur 1: flowschema van dataverwerking in R._


#### Verkregen data

De RNA-sequencing data zijn verkregen uit synoviumbiopten van 4 patiënten met RA en 4 gezonde controles. Bij de RA-patiënten is de diagnose bevestigd door de aanwezigheid van anti-CCP autoantistoffen. Een overzicht van de data is te vinden in de ruwe data link.

_Tabel 1. overzicht van de gebruikte samplemetadata in deze studie._

| Monster ID            | Leeftijd | Geslacht    | Conditie                             |
|-----------------------|----------|-------------|--------------------------------------|
| SRR4785819            | 31       | Vrouw       | Normaal                              |
| SRR4785820            | 15       | Vrouw       | Normaal                              |
| SRR4785828            | 31       | Vrouw       | Normaal                              |
| SRR4785831            | 42       | Vrouw       | Normaal                              |
| *Gemiddelde leeftijd* | *29.8*   | –           | *Normaal*                            |
| SRR4785979            | 54       | Vrouw       | Rheumatoid arthritis (Vastgesteld)   |
| SRR4785980            | 66       | Vrouw       | Rheumatoid arthritis (Vastgesteld)   |
| SRR4785986            | 60       | Vrouw       | Rheumatoid arthritis (Vastgesteld)   |
| SRR4785988            | 59       | Vrouw       | Rheumatoid arthritis (Vastgesteld)   |
| *Gemiddelde Leeftijd* | *59.8*   | –           | *Rheumatoid arthritis (Vastgesteld)* |


#### Mappen van data en countmatrix 
De ruwe sequencingbestanden (FASTQ-formaat) zijn opgeslagen in de map [Ruwe Data](Ruwe%20Data). Voor de uitlijning is het humane referentiegenoom GRCh38 (GCF_000001405.40, versie GRCh38.p14) van NCBI gebruikt. Vanwege de grootte van het genoom is het FASTA-bestand niet opgenomen in de repository, maar kan via NCBI worden gedownload via accessionnummer [GCF_000001405.40](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/).

De referentie-index is opgebouwd met het R-pakket Rsubread (versie 2.20.0) met voldoende geheugen (4 GB) om de indexbestanden te splitsen. De sequencing reads zijn uitgelijnd met de functie align(), waarna de resulterende BAM-bestanden gesorteerd en geïndexeerd zijn met Rsamtools (Versie 2.14.0). De genexpressie is geteld met featureCounts(), waarbij gebruik is gemaakt van de GTF-annotatie Homo_sapiens.GRCh38.114.chr.gtf.gz van ENSEMBL.

Differentiële expressieanalyse is uitgevoerd met DESeq2 (Versie 1.46.0), waarbij het model ~ treatment (controle vs. RA) is toegepast. Resultaten zijn gevisualiseerd met EnhancedVolcano (Versie 1.18.0), en geanalyseerd met goseq (Versie 1.52.0) en pathview (Versie 1.38.0) om biologische functies en pathways te identificeren.

#### Statistische analyse
De gegenereerde countmatrix (data/count_matrix.txt) en een behandel-tabel met controle- en RA-status zijn ingeladen in DESeq2 (versie 1.46.0) om differentiële genexpressie te berekenen. Hierbij zijn log₂ fold changes, p-waarden en meervoudige testcorrecties (Benjamini-Hochberg) bepaald. Resultaten zijn gevisualiseerd in een volcanoplot (log₂ fold change vs. significantie). Daarnaast is functionele verrijkingsanalyse uitgevoerd met goseq voor Gene Ontology en met pathview voor KEGG pathway-analyse.

## Resultaten

### Differentiële genexpressieanalyse 
Om genen te identificeren die verschillend tot expressie kwamen tussen RA-patiënten en gezonde controles, werd een RNA-seq analyse uitgevoerd met DESeq2. In totaal werden 102 genen met verhoogde expressie (log₂FC > 1, padj < 0.05) en 88 genen met verlaagde expressie (log₂FC < -1, padj < 0.05) geïdentificeerd. De resultaten worden weergegeven in een [volcanoPlot](Resultaten/Plots/VolcanoplotRA.png) , waarin de log2 fold change tegenover de -log10 p-waarde van elk gen staat. De volledige lijst van significante genen is beschikbaar in [ResultatenRA](Resultaten/ResultatenRA.csv).

De plot toonde aan dat er meerdere genen significant op- of neerwaarts gereguleerd zijn. Echter, werd in de visualisatie de genen weergegeven als punten, maar zonder Labels (namen van de genen) erbij. Hierdoor was het lastig om direct af te lezen welke specifieke genen betrokken waren bij RA.Voor toekomstige analyse zou het toevoegen van genlabels helpen om de interpretatie te verbeteren. 
Door technische problemen in R kon dit op dit moment niet worden gerealiseerd. 
 
### GO-enrichmentanalyse
Om de biologische betekenis van de differentieel tot expressie gebrachte genen te achterhalen, werd een [GO-enrichmentanalyse](Resultaten/Plots/Rplot03.png) uitgevoerd. De resultaten werden weergegeven in een barplot waarin termen zoals: _Immune system process_ (processen van het immuunsysteem)_, Immune response_ (reactie van het immuunsysteem op prikkels) _, Protein binding_ ( binding van eiwitten aan andere moleculen) _, Intracellular organelle lumen_ (de binnenruimte van celorganellen) duidelijk naar voren komen. Deze termen waren statistisch verrijkt en wezen op verstoring van immuunprocessen, wat goed paste bij het ziektebeeld van RA. 
De visualisatie toonde zowel het aantal betrokken genen als p-waarden, en gaf een helder overzicht van de belangrijkste biologische processen.

### KEGG-pathwayanalyse
Voor verdere interpretatie werd een [KEGG-pathwayanalyse](Resultaten/hsa03260.png) uitgevoerd. Het bijgevoegde figuur toonde een HIV-gerelateerde pathway, wat niet direct relevant leek voor RA. Toch bevatte deze pathway immuunreceptoren zoals _CD4_, die ook een rol spelen in auto-immuunziekten zoals RA. 

Hoewel deze pathway enkele relevante elementen bevat, zou het beter zijn om pathways te tonen die direct met RA te maken hebben, zoals: 
_Cytokine-cytokine receptor interaction, Toll-like receptor signaling pathway_ , deze zijn nauw betrokken bij ontstekingsreacties en geven een duidelijker beeld van de moleculaire processen in RA.


## Conclusie
Deze RNA-sequencinganalyse van synoviumbiopten van RA-patiënten en gezonde controles laat zien dat er duidelijke verschillen zijn in genexpressie tussen beide groepen. De differentiële expressieanalyse identificeerde meerdere genen met significante op- of neerregulatie, hoewel verdere interpretatie beperkt werd door het ontbreken van genlabels in de visualisatie.

De GO-enrichmentanalyse toonde een duidelijke verrijking van immuun-gerelateerde processen, zoals _immune response_ en _Protein binding_, wat goed aansluit bij het ontstekingskarakter van RA. Dit bevestigt dat immuunactivatie een centrale rol speelt in de pathogenese van de ziekte.

Hoewel de KEGG-pathwayanalyse oorspronkelijk een HIV-gerelateerde pathway toonde, bevatte deze wel immuunreceptoren die ook relevant zijn voor RA. voor een meer directe interpretatie zouden RA-specifieke pathways zoals _Cytokine-cytokine receptor interaction_ beter geschikt zijn. 

Afsluitend bieden de resultaten waardevolle inzichten in de moleculaire processen die betrokken zijn bij RA, met name op het gebied van immuunregulatie. Verdere verfijning van de visualisaties en pathwayselectie kan bijdragen aan een nog scherpere interpretatie en mogelijke aanknopingspunten voor therapie of diagnostiek.




