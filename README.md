

## Inhoud/structuur

- `Bam bestanden/` - reads dat gemapt zijn tegen referentiegenoom  
- `bronnen/` - gebruikte bronnen -> note: toevoegen NCBI in wordt document
- `RScripts/` – Rscript waar alle data wordt verwerkt 
- `resultaten/` - grafieken en tabellen
- `Ruwe Data/` – Hierin zit de data van de Patiënten met en zonder Reumatoïde Atritis (RA).
- `assets/` - overige documenten voor de opmaak van deze pagina
- `README.md` - Is het document waar het verslag (de tekst) in wordt gemaakt. 






---

## Inleiding

Reumatoïde artritis (RA) is een chronische, systematische auto-immuunziekte waarbij het immuunsysteem het eigen lichaam aanvalt. De exacte oorzaak is nog onbekend, maar vermoedelijk speelt een combinatie van genetische aanleg, omgevingsfactoren en een ontspoord immuunsysteem een rol [Gabriel, 2001](https://doi.org/10.1016/S0889-857X(05)70201-5). Een belangrijk kenmerk van RA is synovitis: een ontsteking van het gewrichtsslijmvlies, wat leidt tot pijn, zwelling en uiteindelijk gewrichtsschade [Radu & Bungau, 2021](https://doi.org/10.3390/cells10112857). Vroege opsporing en behandeling kunnen schade beperken, maar genezing is (nog) niet mogelijk.

In dit project wordt RNA-sequencing data geanalyseerd van synoviumbiopten van zowel gezonde personen als patiënten met RA (gevestigde diagnose, > 12 maanden). De analyse is uitgevoerd in R, met als doel: inzicht krijgen in genexpressieverschillen tussen de groepen en achterhalen welke biologische processen betrokken zijn bij RA, via Gene Ontology (GO-analyse) 

De gebruikte brondata en artikelen zijn te vinden in de folder [bronnen](Bronnen). 



## Methoden

Voor deze studie is RNA-sequencing data geanalyseerd van synoviumbiopten van gezonde controles en patiënten met reumatoïde artritis (RA).

*Figuur 1: flowschema van dataverwerking in R. nog maken*


#### Verkregen data

De RNA-sequencing data zijn verkregen uit synoviumbiopten van 4 patiënten met RA en 4 gezonde controles. Bij de RA-patiënten is de diagnose bevestigd door de aanwezigheid van anti-CCP autoantistoffen. Een overzicht van de data is te vinden in  ruwe data link(data link ruwe data).

#### Mappen van data en countmatrix 
De ruwe sequencingbestanden (FASTQ-formaat) zijn opgeslagen in de map Ruwe Data/. Voor de uitlijning is het humane referentiegenoom GRCh38 (GCF_000001405.40, versie GRCh38.p14) van NCBI gebruikt. Vanwege de grootte van het genoom is het FASTA-bestand niet opgenomen in de repository, maar kan via NCBI worden gedownload via accessionnummer GCF_000001405.40.

De referentie-index is opgebouwd met het R-pakket Rsubread (versie 2.20.0) met voldoende geheugen (4 GB) om de indexbestanden te splitsen. De sequencing reads zijn uitgelijnd met de functie align(), waarna de resulterende BAM-bestanden gesorteerd en geïndexeerd zijn met Rsamtools. De genexpressie is geteld met featureCounts(), waarbij gebruik is gemaakt van de GTF-annotatie Homo_sapiens.GRCh38.114.chr.gtf.gz van ENSEMBL.

Differentiële expressieanalyse is uitgevoerd met DESeq2, waarbij het model ~ treatment (controle vs. RA) is toegepast. Resultaten zijn gevisualiseerd met EnhancedVolcano en geanalyseerd met goseq en pathview om biologische functies en pathways te identificeren.

#### Statistische analyse
De gegenereerde countmatrix (data/count_matrix.txt) en een behandel-tabel met controle- en RA-status zijn ingeladen in DESeq2 (versie 1.46.0) om differentiële genexpressie te berekenen. Hierbij zijn log₂ fold changes, p-waarden en meervoudige testcorrecties (Benjamini-Hochberg) bepaald. Resultaten zijn gevisualiseerd in een volcano plot (log₂ fold change vs. significatie). Daarnaast is functionele verrijkingsanalyse uitgevoerd met goseq voor Gene Ontology en met pathview voor KEGG pathway-analyse.

## Resultaten
Resultaten: +- 200 woorden, inclusief correcte verwijzingen.
Resultaten die verkregen zijn uit R

Voorbeeld:

De analyse van de RNA-seq data toonde duidelijke verschillen in genexpressie tussen reumatoïde artritis (RA)-patiënten en gezonde controles. In totaal werden 102 genen met verhoogde expressie (log₂FC > 1, padj < 0.05) en 88 genen met verlaagde expressie (log₂FC < -1, padj < 0.05) geïdentificeerd. De volledige lijst van significante genen is beschikbaar in het resultatenbestand (welke okalweer) 


De [volcano](Resultaten/Plots/VolcanoplotRA.png) plot toont een duidelijke scheiding tussen genen met verhoogde en verlaagde expressie in RA-patiënten. Genen met de meest significante expressieveranderingen (padj < 0.05) vallen op door hun sterke betrokkenheid bij ontstekingsgerelateerde processen.

Vervolgonderzoek met GO-enrichmentanalyse toonde een significante oververtegenwoordiging van biologische processen die samenhangen met het immuunsysteem, waaronder ‘T cel activatie’ en ‘cytokine-mediated signaling’ (go resultaten LINKJE). KEGG-pathwayanalyse(linkje) suggereerde dat vooral routes zoals ‘T cell receptor signaling’ en ‘cytokine-cytokine receptor interaction’ betrokken zijn bij RA.
Alle benodigde data, R-scripts en visualisaties zijn georganiseerd en terug te vinden in de betreffende mappen van de GitHub-repository 


## Conclusie
Conclusie: +- 200 woorden, inclusief aanbevelingen en onderzoek in context
plaatsen.

In dit onderzoek is RNA-sequencing gebruikt om genexpressieverschillen te analyseren tussen synoviumbiopten van gezonde controles en patiënten met reumatoïde artritis (RA). De differentiële expressieanalyse toonde aan dat een aanzienlijk aantal genen significant anders tot expressie komt in RA, waarbij met name genen betrokken bij immuunrespons en ontstekingsprocessen zijn op- of afgeschaald. Dit bevestigt het ontstekingskarakter van RA zoals beschreven in de literatuur en sluit aan bij het klinische beeld van synovitis, waarbij het gewrichtsslijmvlies ontstoken raakt en beschadigd wordt.

De verrijkingsanalyse via Gene Ontology en KEGG benadrukte het belang van pathways die gerelateerd zijn aan immuunactivatie, cytokinesignalisatie en celadhesie, wat het mechanisme van de ontstekingsreactie bij RA verder ondersteunt. Deze bevindingen dragen bij aan het inzicht in de moleculaire processen die bijdragen aan de progressie van RA en kunnen potentieel leiden tot nieuwe biomarkers voor vroege diagnose of doelgerichte therapieën.

Voor toekomstig onderzoek wordt aanbevolen om deze genen en pathways verder te valideren met experimentele technieken en om grotere patiëntengroepen te analyseren om de klinische relevantie te bevestigen




