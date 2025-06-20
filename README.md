# Transcriptomics analyse reumatoide-artritis
Transcriptomics analyse in R van eiwit-gerelateerde genexpressie in reumatoïde artritis versus gezonde controles 

<figure style="text-align: center;">
  <img src="afbeelding_voorblad/normal_vs_reuma.png" alt="Afbeelding voorblad" width="800" style="margin: auto; display: block;">
</figure>

[(Ding. et al., 2023)](bronnen/ding_etal_2023.pdf)

# Introductie
Reumatoïde artritis (RA) is een chronische auto-immuun ontstekingsziekte en brengt schade en functieverlies toe aan de gewrichten. De ziekte kan zich ook extra-articulair manifesteren, waarbij er aantasting op treed aan de huid, ogen, hart, longen, nieren, zenuwstelsel en maagdarmstelsel [(Radu & Bungau, 2021)](bronnen/radu_etal_2021.pdf). RA komt wereldwijd voor bij 1-2% van de bevolking, en wordt bij vrouwen twee tot drie keer vaker aan getroffen dan in mannen [(Kang et al., 2022)](bronnen/kang_etal_2022.pdf). RA ontwikkeld zich in het begin met aspecifieke symptomen die kunnen overlappen met andere ziekten, dit maakt het lastig de aandoening in een vroeg stadium vast te stellen. Het niet behandelen van de ziekte leidt tot een verhoogd functieverlies en mortaliteit [(Chauhan., et al)](bronnen/chauchan_etal_2023.pdf). Verschillende medicijnen worden toegepast om de ziekte onder controle te houden en richten zich op ontstekingsmediatoren zoals tumornecrosefactor (TNF)-α, interleukine (IL)-6, en B-cellen. Deze behandelingen zijn niet voor iedere RA patiënt effectief, maar leiden wel tot nadelige bijwerkingen [(Wang et al., 2024)](bronnen/wang_etal_2024.pdf). Er is meer onderzoek nodig naar het identificeren van RA in patiënten en betere behandeling met geneesmiddelen. In dit onderzoek wordt een transcriptomics analyse uitgevoerd op RNA- sequence gegevens van RA patiënten. Hiermee wordt het genexpressie profiel van RA patiënten in kaart gebracht om
biomarkers optesporen, en om mogelijk bij te dragen aan een beter begrip van potentiële behandelmogelijkheden, cel regulatie en de regulerende netwerken in de ontwikkeling van RA.


# Methode
**Synoviale weefselmonsters**

De RNA-sequence gegevens van de synovium biopten werden verkregen uit eerder uitgevoerd onderzoek [(Platzer et al., 2019)](bronnen/platzer_etal_2019.pdf). In totaal werden 8 monsters verzameld, bestaande uit monsters van 4 ACPA positieve vrouwen met RA (gemiddelde leeftijd 59.8 ± 4.9) en monsters van 4 gezonde vrouwen (gemiddelde leeftijd 29.8 ± 11.1 ). 

**Data analyse**

Data analyse van de gegevens werd uitgevoerd in R studio. Reads zijn aan de hand van kwaliteit controles getrimd en vervolgens gemapt tegen het humane referentiegenoom Homo_sapiens.GRCh38.114 (ensemble) met behulp van de align() functie uit de `Rsubread (V2.22.1)` package, te zien in [script_mappen](scripts/mapping_day_1.R). Genexpressie-kwantificatie werd uitgevoerd met featureCounts functie uit de Rsubread package, resulterend in een gen-telling matrix. Significant verschillen in genexpressie werden geanalyseerd met de `DESeq2 (V1.48.1)` package, te zien in [script_DESeqanalyse](scripts/Analyse_en_statistiek_day_3.R). Significante genen uit DESeq analyse werden gevisualiseerd met een volcanoplot, met behulp van de package `EnhancedVolcano (V1.26.0)`.  Om het aantal fout-positieve resultaten te beperken werden de p-waardes (FDR), gecorrigeerd met de Benjamini-Hochberg methode. KEGG enrichment analyse werd uitgevoerd met de DESeq resultaten om de meest verijkte pathway te identificeren. Om inzicht te krijgen in de genen van de meest verijkte pathway werd er een KEGG-analyse uitgevoerd met de package `KEGGREST (V1.48.0)`, te zien in [script_KEGG_analyse](scripts/Analyse_en_statistiek_day_3.R). Gene Ontology (GO)-enrichment analyse werd uitgevoerd om inzicht te krijgen in de meest verijkte biologische processen, met de packages `clusterProfiler (V4.16.0)`, `org.Hs.eg.db (V3.21.0)`, `enrichplot (V1.28.1)`. Onderscheid werd gemaakt tussen verlaagde en verhoogde genexpressie met een log2 fold change drempel van 0.5 en -0.5. De top 5 verhoogde en verlaagde GO-termen werden gevisualiseerd met de package `ggplot2 (V3.5.2)`, te zien in [script_GO_analyse](scripts/GO_analysis_script.R). 

Toevoegen plaatje workflow

# Resultaten 


**DESeq analyse**

Significante genen uit de DESeq analyse werden gevisualiseerd, te zien in de volcanoplot. Log2 fold change op de x- as uitgezet tegen -Log₁₀ P op de y-as. De beteknis van de kleuren zijn: grijs, niet significant,  groen, alleen de fold change is significant, rood, zowel de verandering in expressie als de p-waarde zijn significant. De rood gekleurde genen uit de volcano plot en hun mogelijke functie staan in tabel 1.

*Tabel 1; top 5 meest significante genen uit de volcano plot. Gebasseerd op een hoge count, een hoge fold enrichment en een lage p.adjust en p-value. Functie van genen wordt weergeven.* 

| Rang | Gen       | Functie / Betekenis |
|------|-----------|---------------------|
| 1    | ANKRD30BL | Mogelijk betrokken bij signaaltransductie en celprocessen (Almeida. et al., 2020) |
| 2    | MT-ND6    | Onderdeel van mitochondriaal complex I, essentieel voor energieproductie. |
| 3    | ZNF598    | Speelt een rol in kwaliteitscontrole van ribosomen. |
| 4    | CROCC     | Belangrijk voor cohesie van het centrosoom tijdens celdeling. |
| 5    | IKBKG     | Betrokken bij regulatie van immuunresponsen. |

De pathways waarin de meeste genen uit DESeq analyse voorkwamen zijn te zien in de dotplot. Op de y-as de verrijkte biologische pathways en op de y-as de genratio. De top 3 pathways waren de MAPK signaling pathway, TNF signaling pathway en NF-kappa B signaling pathway. NF-kappa B signaling pathway werd als een van de belangrijkste pathways geselecteerd, omdat de pathway consistent hoog scoort op meerdere statistische criteria, te zien in de KEGG enrichment analyse. De criteria bestonden uit een hoge count, een hoge fold enrichment en een lage p.adjust en p-value.

**NF-kappa B signaling pathway**

De genen uit de NF-kappa B signaling pathway werden gevisualiseerd in een signaalroutekaart, te zien in de kaart. De genen of eiwitten met verhoogde expressie werden in rood weergegeven en de genen of eiwitten met verlaagde expressie in het groen. De kleurenschaal loopt van -1 (groen) tot +1 (rood). Hierin was te zien dat er een verhoogde expressie van de genen RELA (p65), NFKB1 (p50) en IKBKG (NEMO) aanwezig was in de  NF-kappa B signaling pathway. Dit komt overeen met wat bekend is in de literatuur waarbij een hoge expressie van deze genen worden aangetroffen in RA patienten (Nejatbakhsh. et al., 2020).  

**Gene Ontology**

De top 5 verrijkte GO-termen voor verhoogde en verlaagde genexpressie zijn weergegeven in de [barplot](resultaten/Top5_GOtermen.png). GO term *immune response-regulating cell surface receptor signaling pathway* bleek het meest verijkt te zijn met een verhoogde expressie van genen. Dit komt overeen met wat bekend is uit de litertuur, waarin CD4+ T-cellen, pathogene B-cellen, macrofagen, ontstekingscytokinen, chemokinen en autoantilichamen autoreactief reageren bij RA patienten (Jang. et al., 2022). Verlaagde expressie van genen werdt voornamelijk aangetroffen in de GO-term  *pattern specification process*. In de literatuur is ook terug te lezen dat RA invloed heeft op de processen voor het vormen van hyperplastisch synovium, kraakbeenschade en boterosie. Deze processen zijn afhankelijk van signaalroutes die ook betrokken zijn bij de GO-term (Guo. et al., 2018).

# conclusie 
Benomen genen KEGG pathway
Benoemen meest verijkt genen in GO-term

Voor de resultaten van je eigen onderzoek, die je in het verleden hebt uitgevoerd, gebruik je de onvoltooid verleden tijd. 
Voor algemene bevindingen of bekende feiten, kun je de onvoltooid tegenwoordige tijd gebruiken

