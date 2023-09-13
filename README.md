## Search for Z-DNA regions and quadruplexes in gene promoters of representatives of the genus *Plasmodium* (type *Apicomplexans*)

[This README in russian](README.ru.md)

This work is part of a larger analysis presented [here](https://github.com/OMALOE/hse22_project_Apicomplexans).

Z-DNA is left-handed DNA in a specific conformation. Its functions are not fully learnt, but among [other things](https://www.imrpress.com/journal/FBL/12/7/10.2741/2399), Z-DNA may be involved in the regulation of gene expression. This work represents an attempt to search for proteins whose expression could potentially be regulated by Z-DNA.

For this analysis were taken 5 species that cause malaria:
* *P. falciparum* - in humans
* *P. vivax* - in humans
* *P. knowlesi* - in humans and some macaques
* *P. gaboni* - in chimpanzees
* *P. yoelii* - in rodents

Data for these five organisms were downloaded from [here](https://www.ncbi.nlm.nih.gov/genome/browse#!/eukaryotes/refseq_category:representative) and [here](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA).  
All .fasta and .gb (INSDC) files were downloaded from the NCBI browser using a bash script. Total 70 .fasta files and 70 .gb files (one file per chromosome).  
At the same time files feature_table.txt and protein.faa (for compiling protein clusters) were downloaded via ftp.

In addition to the python code (src folder), the Linux terminal was actively used (some commands in the .docx file).

So, the regions with Z-DNA were detected using the Z-Hunt algorithm (the longer sequence, the higher ZH-Score).

|Species|Number of chromosomes|Genome length|GC-content|Number of annotated genes|Fraction of annotated genes|Number of regions with significant ZH_SCORE*|Length of sections with significant ZH-SCORE*|
|-|-|-|-|-|-|-|-|
|*P. falciparum*|14**|23 292 622|19.34|5387|0.59|1064|13 714|
|*P. vivax*|14|22 621 071|44.87|5389|0.62|37018|476 891|
|*P. gaboni*|14|18 878 758|17.84|5198|0.65|617|8 162|
|*P. knowlesi*|14**|23 938 832|38.61|5326|0.59|18210|234 759|
|*P. yoelii*|14**|23 043 114|21.74|6039|0.56|1844|24 712|

`*` Namely, with ZH-SCORE >= 500.  
`**` In addition to 14 chromosomes, there were also data on the mitochondrial genome and apicoplast genome, but there was only one region with ZH-SCORE >= 500, so they were not considered.  

Using intersectBed, Z-DNA regions were intersected with 600 nucleotide ranges around the start of the ORF (300 nucleotides in each direction).  
If Z-DNA region fell there, then it was believed that a relationship with its regulation may be possible. If this region was before the start of the ORF, then it was roughly called a hit in the promoter. If after, it was roughly called a hit in the gene (without division into exons/introns).

To search for conservation within the genus, protein clustering was applied.  
This was done using proteinortho (based on the BLAST algorithm).

#### For the species above, proteinortho constructed 5319 clusters.

The number of clusters consisting of n genomes had the following distribution (statistics on the number of orthologs in a cluster):

![clust_by_number_of_genoms](https://user-images.githubusercontent.com/60808642/173706235-3e63ee97-388a-4c93-91a8-9709d930ce60.png)

The number of clusters, which include n protein sequences, had the following distribution (sometimes there were several paralogues for one species):

![clust_by_number_of_proteins](https://user-images.githubusercontent.com/60808642/173706248-663fa365-f091-4f64-818e-d7a863e28d49.png)

I selected 10 clusters for closer examination. They were selected first of all by the presence of Z-DNA near the beginning of the gene for most species in a cluster, and then by the total ZH-SCORE.  
Unfortunately, for only two clusters there was Z-DNA near the ORF for all five species. The remaining clusters shown here contain Z-DNA for 4 of the 5 species.

For each cluster, a protein alignment was also constructed (using the ClustalW algorithm in the MEGAX program). All files related to this can be found in the clusters folder; below are screenshots of the beginning and end of the protein alignment (for cluster 9, the protein turned out to be short enough that it fit entirely into one screenshot).

#### Cluster 1
|Species|Proteins in the cluster|IDs of encoded proteins|Functions of encoded proteins|Z-DNA location|ZH_SCORE|
|-|-|-|-|-|-|
|*P. falciparum*|1|CZT99390.1|cAMP-dependent protein kinase regulatory subunit|Gene|583.4285|
|*P. vivax*|1|EDL47592.1|cAMP-dependent protein kinase regulatory subunit, putative|Promoter|138924.1|
|*P. gaboni*|1|KYN98034.1|cAMP-dependent protein kinase regulatory subunit|Gene|583.4285|
|*P. knowlesi*|1|CAA9991018.1|cAMP-dependent protein kinase regulatory subunit, putative|Gene|583.4285|
|*P. yoelii*|1|VTZ81713.1|cAMP-dependent protein kinase regulatory subunit, putative|Gene|583.4285|

![1_start](https://user-images.githubusercontent.com/60808642/174448980-2f576943-8b34-4642-bbf4-f8871227cb43.png)
![1_end](https://user-images.githubusercontent.com/60808642/174448983-2682eaed-1240-4730-aa39-0612a8875a66.png)
![cluster1](https://user-images.githubusercontent.com/60808642/174160532-0c6192f2-66b1-4b3b-a79f-9f6f77b8811c.png)

#### Cluster 2
|Species|Proteins in the cluster|IDs of encoded proteins|Functions of encoded proteins|Z-DNA location|ZH_SCORE|
|-|-|-|-|-|-|
|*P. falciparum*|1|CAD51935.1|protein phosphatase-beta|No Z-DNA|0.0|
|*P. vivax*|1|EDL45060.1|protein phosphatase-beta, putative|Promoter (outside chart)|94590.41|
|*P. gaboni*|1|KYO00236.1|protein phosphatase-beta|Promoter|705.4245|
|*P. knowlesi*|1|CAA9987523.1|protein phosphatase-beta, putative|Promoter|650.9198|
|*P. yoelii*|1|VTZ77586.1|protein phosphatase-beta, putative|Promoter|583.4285|

![2_start](https://user-images.githubusercontent.com/60808642/174449009-afaf75c1-9997-445e-b20c-3b6465ce51dc.png)
![2_end](https://user-images.githubusercontent.com/60808642/174449010-c836b79b-61be-45d3-a472-89ab93af49b6.png)
![cluster2](https://user-images.githubusercontent.com/60808642/174160560-c4cf8210-00d5-47ae-ae28-40608a534ed4.png)

#### Cluster 3
|Species|Proteins in the cluster|IDs of encoded proteins|Functions of encoded proteins|Z-DNA location|ZH_SCORE|
|-|-|-|-|-|-|
|*P. falciparum*|1|CAB11145.2|formate-nitrite transporter|Gene|2091.083|
|*P. vivax*|1|EDL44820.1|transporter, putative|Promoter (outside chart)|38833.58|
|*P. gaboni*|1|KYO03279.1|putative formate-nitrite transporter|Gene|6565.992|
|*P. knowlesi*|1|CAA9987895.1|formate-nitrite transporter, putative|Gene|28780.5|
|*P. yoelii*|1|VTZ73436.1|formate-nitrite transporter, putative|No Z-DNA|0.0|

![3_start](https://user-images.githubusercontent.com/60808642/174449018-7a893939-b6d4-4e20-b9a6-2f787a5dd5bf.png)
![3_end](https://user-images.githubusercontent.com/60808642/174449021-fe072a9e-153b-4542-82c8-41c20e011986.png)
![cluster3](https://user-images.githubusercontent.com/60808642/174160577-2a079fbd-54cf-4932-a4d4-49fbbf84ca5a.png)

#### Cluster 4
|Species|Proteins in the cluster|IDs of encoded proteins|Functions of encoded proteins|Z-DNA location|ZH_SCORE|
|-|-|-|-|-|-|
|*P. falciparum*|1|CZT98040.1|octaprenyl pyrophosphate synthase|Promoter|1743.107|
|*P. vivax*|1|EDL43227.1|prenyl transferase, putative|Promoter|16770.32|
|*P. gaboni*|1|KYO03461.1|octaprenyl pyrophosphate synthase|No Z-DNA|0.0|
|*P. knowlesi*|1|CAA9986822.1|octaprenyl pyrophosphate synthase, putative|Promoter|2180.071|
|*P. yoelii*|1|VTZ72360.1|octaprenyl pyrophosphate synthase, putative|Promoter|583.4285|

![4_start](https://user-images.githubusercontent.com/60808642/174449028-b40aa99f-8bfd-464a-b99d-93556b9e6fa5.png)
![4_end](https://user-images.githubusercontent.com/60808642/174449030-9036d647-9bed-409d-91b3-ce5604ca6be3.png)
![cluster4](https://user-images.githubusercontent.com/60808642/174160605-0062db79-4ad6-49ff-9265-e3f7d19497b8.png)

#### Cluster 5
|Species|Proteins in the cluster|IDs of encoded proteins|Functions of encoded proteins|Z-DNA location|ZH_SCORE|
|-|-|-|-|-|-|
|*P. falciparum*|1|CAD52619.1|peptidase, putative|Gene|612.3848|
|*P. vivax*|1|EDL44548.1|hypothetical protein, conserved|Gene|8323.257|
|*P. gaboni*|1|KYN97268.1|putative membrane protein|Gene|612.3848|
|*P. knowlesi*|1|CAA9989837.1|peptidase, putative|Gene|8323.257|
|*P. yoelii*|1|VTZ81050.1|peptidase, putative|No Z-DNA|0.0|

![5_start](https://user-images.githubusercontent.com/60808642/174449039-128108f5-c6a9-4130-be26-8ecd8e7cd876.png)
![5_end](https://user-images.githubusercontent.com/60808642/174449041-fb0e31ee-2144-47bc-82a5-8a908c07a412.png)
![cluster5](https://user-images.githubusercontent.com/60808642/174160640-7beb2aa8-bd5f-4387-9bb1-960a76760d33.png)

#### Cluster 6
|Species|Proteins in the cluster|IDs of encoded proteins|Functions of encoded proteins|Z-DNA location|ZH_SCORE|
|-|-|-|-|-|-|
|*P. falciparum*|1|CZT98598.1|conserved protein, unknown function|Promoter|915.9191|
|*P. vivax*|1|EDL42405.1|hypothetical protein PVX_111020|Promoter|13713.99|
|*P. gaboni*|1|KYN99721.1|hypothetical protein PGSY75_1034100|Promoter|915.9191|
|*P. knowlesi*|1|CAA9987243.1|conserved protein, unknown function|Promoter|908.3955|
|*P. yoelii*|1|VTZ74384.1|conserved protein, unknown function|No Z-DNA|0.0|

![6_start](https://user-images.githubusercontent.com/60808642/174449046-bf327b64-b73e-4d61-bd19-689173c8e725.png)
![6_end](https://user-images.githubusercontent.com/60808642/174449047-d76e78c2-dfb1-4191-8e69-4da8a2332b68.png)
![cluster6](https://user-images.githubusercontent.com/60808642/174160667-bd6e244a-7c81-4557-8677-39a5a4bdf0b4.png)

#### Cluster 7
|Species|Proteins in the cluster|IDs of encoded proteins|Functions of encoded proteins|Z-DNA location|ZH_SCORE|
|-|-|-|-|-|-|
|*P. falciparum*|1|CZT98374.1|hypoxanthine-guanine phosphoribosyltransferase|Gene|583.4285|
|*P. vivax*|1|EDL44708.1|hypoxanthine phosphoribosyltransferase, putative|Promoter|2276.8|
|*P. gaboni*|1|KYN99499.1|hypoxanthine-guanine phosphoribosyltransferase|No Z-DNA|0.0|
|*P. knowlesi*|1|CAA9987765.1|hypoxanthine-guanine phosphoribosyltransferase, putative|Gene|6408.923|
|*P. yoelii*|1|VTZ79941.1|hypoxanthine-guanine phosphoribosyltransferase, putative|Gene|2091.083|

![7_start](https://user-images.githubusercontent.com/60808642/174449051-114e20f3-4b1f-436c-ae4c-8d9c6404e025.png)
![7_end](https://user-images.githubusercontent.com/60808642/174449052-ed5e628c-89ef-4a53-b453-d8f8ae18d23b.png)
![cluster7](https://user-images.githubusercontent.com/60808642/174160691-a5199f69-0862-40c3-9eb1-10568d1bc7f5.png)

#### Cluster 8
|Species|Proteins in the cluster|IDs of encoded proteins|Functions of encoded proteins|Z-DNA location|ZH_SCORE|
|-|-|-|-|-|-|
|*P. falciparum*|1|CZT98084.1|conserved Plasmodium protein, unknown function|Promoter|2183.574|
|*P. vivax*|1|EDL43268.1|hypothetical protein, conserved|Promoter|1820.585|
|*P. gaboni*|1|KYO03505.1|hypothetical protein PGSY75_0207100|Promoter|2183.574|
|*P. knowlesi*|1|CAA9986780.1|conserved Plasmodium protein, unknown function|Promoter|4831.423|
|*P. yoelii*|1|VTZ72434.1|conserved Plasmodium protein, unknown function|No Z-DNA|0.0|

![8_start](https://user-images.githubusercontent.com/60808642/174449060-4b59cfbc-c4f4-4504-90ab-f20253d21ef8.png)
![8_end](https://user-images.githubusercontent.com/60808642/174449063-2d49adc8-1eac-4985-a54f-837b118fa74a.png)
![cluster8](https://user-images.githubusercontent.com/60808642/174160707-d2d6a252-a5bf-48c4-9459-0d02bda79429.png)

#### Cluster 9
|Species|Proteins in the cluster|IDs of encoded proteins|Functions of encoded proteins|Z-DNA location|ZH_SCORE|
|-|-|-|-|-|-|
|*P. falciparum*|1|CZT98918.1|conserved protein, unknown function|Promoter|731.2843|
|*P. vivax*|1|EDL45645.1|hypothetical protein PVX_091900|Gene|3661.11|
|*P. gaboni*|1|KYN93456.1|hypothetical protein PGSY75_0013800|No Z-DNA|0.0|
|*P. knowlesi*|1|CAA9988282.1|conserved protein, unknown function|Gene|3661.11|
|*P. yoelii*|1|VTZ78436.1|conserved protein, unknown function|Gene|2173.083|

![9](https://user-images.githubusercontent.com/60808642/174449073-43846759-d111-4eec-b75d-dc2c6cb64f8a.png)
![cluster9](https://user-images.githubusercontent.com/60808642/174160723-50ba630e-251d-42d7-b499-bf1aba398401.png)

#### Cluster 10
|Species|Proteins in the cluster|IDs of encoded proteins|Functions of encoded proteins|Z-DNA location|ZH_SCORE|
|-|-|-|-|-|-|
|*P. falciparum*|1|CAD51193.2|zinc finger protein, putative|Gene|2183.574|
|*P. vivax*|1|EDL45264.1|D13 protein, putative|Gene|826.8355|
|*P. gaboni*|1|KYO00773.1|putative zinc finger protein|Gene|2183.574|
|*P. knowlesi*|1|CAA9986906.1|zinc finger protein, putative|Gene|826.8355|
|*P. yoelii*|1|VTZ76436.1|zinc finger protein, putative|Gene|2183.574|

![10_start](https://user-images.githubusercontent.com/60808642/174449077-a55724f4-6aba-43f6-a3fc-795a08bbac2c.png)
![10_end](https://user-images.githubusercontent.com/60808642/174449078-f1c5a51a-a474-40ba-a689-92e7857d89be.png)
![cluster10](https://user-images.githubusercontent.com/60808642/174160733-f5616104-11ce-4e09-a292-f43eee7be170.png)

### Search for quadruplexes
To find quadruplexes, the pqsfinder library and R code (in the src folder) were used.  
After obtaining the quadruplex SCORE, these regions were intersected in the same way as Z-DNA.  
Apparently, there were even fewer quadruplexes and they were unevenly distributed.  
Thus, in plasmodium knowlesi and plasmodium vivax there were tens of times more of them (which correlates with the amount of Z-DNA in these two species).

#### Cluster 1
|Species|Proteins in the cluster|IDs of encoded proteins|Functions of encoded proteins|Location of quadruplexes|SCORE|
|-|-|-|-|-|-|
|*P. falciparum*|1|CZT99927.1|mitochondrial ribosomal protein L21 precursor, putative|No quadruplexes|0.0|
|*P. vivax*|1|EDL47013.1|hypothetical protein, conserved|Promoter|61.0|
|*P. gaboni*|1|KYN95826.1|putative mitochondrial ribosomal protein L21 precursor|Gene|60.0|
|*P. knowlesi*|1|CAA9990390.1|ribosomal protein L21, mitochondrial, putative|Gene|60.0|
|*P. yoelii*|1|VTZ78962.1|mitochondrial ribosomal protein L21 precursor, putative|No quadruplexes|0.0|

![1_start](https://user-images.githubusercontent.com/60808642/174450423-56ebb836-7a8e-4f9c-be99-139611bb092d.png)
![1_end](https://user-images.githubusercontent.com/60808642/174450425-8c94a2a5-554a-491c-8b3f-df4eb6d553b8.png)
![cluster1](https://user-images.githubusercontent.com/60808642/174450046-05419667-2a3d-4226-9b1b-c6990a96eb14.png)

#### Cluster 2
|Species|Proteins in the cluster|IDs of encoded proteins|Functions of encoded proteins|Location of quadruplexes|SCORE|
|-|-|-|-|-|-|
|*P. falciparum*|1|CAD52695.1|DNA helicase, putative|No quadruplexes|0.0|
|*P. vivax*|1|EDL46632.1|helicase, putative|Gene|97.0|
|*P. gaboni*|1|KYN97338.1|putative DNA helicase|No quadruplexes|0.0|
|*P. knowlesi*|1|CAA9988969.1|DNA helicase, putative|Promoter|60.0|
|*P. yoelii*|1|VTZ79625.1|DNA helicase, putative|No quadruplexes|0.0|

![2_start](https://user-images.githubusercontent.com/60808642/174450433-c779e63b-22c1-47c8-81a9-b19a9fd10a52.png)
![2_end](https://user-images.githubusercontent.com/60808642/174450435-14c43814-4115-4a78-b4fc-4d06248f26af.png)
![cluster2](https://user-images.githubusercontent.com/60808642/174450073-dd3d67aa-76c3-479a-802b-24753a731052.png)

#### Cluster 3
|Species|Proteins in the cluster|IDs of encoded proteins|Functions of encoded proteins|Location of quadruplexes|SCORE|
|-|-|-|-|-|-|
|*P. falciparum*|1|CZT98927.1|WD repeat-containing protein, putative|No quadruplexes|0.0|
|*P. vivax*|1|EDL45653.1|hypothetical protein, conserved|Promoter (outside chart)|73.0|
|*P. gaboni*|1|KYN98904.1|hypothetical protein PGSY75_1126500|No quadruplexes|0.0|
|*P. knowlesi*|1|CAA9988291.1|WD repeat-containing protein, putative|Promoter|73.0|
|*P. yoelii*|1|VTZ78428.1|WD repeat-containing protein, putative|No quadruplexes|0.0|

![3_start](https://user-images.githubusercontent.com/60808642/174450440-02e1dbbb-8006-4abc-bcf2-71a0fe97d35c.png)
![3_end](https://user-images.githubusercontent.com/60808642/174450441-3009365a-a644-4011-a2e4-3a23ca591f32.png)
![cluster3](https://user-images.githubusercontent.com/60808642/174450079-1fe90126-c799-460b-9215-269983640ab3.png)

#### Cluster 4
|Species|Proteins in the cluster|IDs of encoded proteins|Functions of encoded proteins|Location of quadruplexes|SCORE|
|-|-|-|-|-|-|
|*P. falciparum*|1|CZT99882.1|conserved Plasmodium protein, unknown function|No quadruplexes|0.0|
|*P. vivax*|1|EDL47055.1|hypothetical protein, conserved|Promoter|64.0|
|*P. gaboni*|1|KYN95782.1|hypothetical protein PGSY75_1417100|No quadruplexes|0.0|
|*P. knowlesi*|1|CAA9990435.1|conserved Plasmodium protein, unknown function|Promoter|73.0|
|*P. yoelii*|1|VTZ79005.1|conserved Plasmodium protein, unknown functio|No quadruplexes|0.0|

![4_start](https://user-images.githubusercontent.com/60808642/174450444-a9e6371f-97cb-4614-aaa5-3b782afe95fd.png)
![4_end](https://user-images.githubusercontent.com/60808642/174450448-9ffb8478-eee9-4819-bf3a-4c873db73a3f.png)
![cluster4](https://user-images.githubusercontent.com/60808642/174450084-d6b0a90e-8582-489b-b68d-42f3fea679b2.png)

#### Cluster 5
|Species|Proteins in the cluster|IDs of encoded proteins|Functions of encoded proteins|Location of quadruplexes|SCORE|
|-|-|-|-|-|-|
|*P. falciparum*|1|CZT99834.1|conserved Plasmodium protein, unknown function|No quadruplexes|0.0|
|*P. vivax*|1|EDL47103.1|hypothetical protein, conserved|Gene|63.0|
|*P. gaboni*|1|KYN95737.1|hypothetical protein PGSY75_1412400|No quadruplexes|0.0|
|*P. knowlesi*|1|CAA9990484.1|conserved Plasmodium protein, unknown function|Gene|73.0|
|*P. yoelii*|1|VTZ79051.1|conserved Plasmodium protein, unknown function|No quadruplexes|0.0|

![5_start](https://user-images.githubusercontent.com/60808642/174450452-1297cb8f-f3df-430e-81bf-a568ec0320c3.png)
![5_end](https://user-images.githubusercontent.com/60808642/174450453-36a8745b-37d0-40d2-93de-1886e17319b3.png)
![cluster5](https://user-images.githubusercontent.com/60808642/174450086-d2a98e6c-9103-473e-9788-8231c96a48d5.png)

#### Cluster 6
|Species|Proteins in the cluster|IDs of encoded proteins|Functions of encoded proteins|Location of quadruplexes|SCORE|
|-|-|-|-|-|-|
|*P. falciparum*|1|CAG25249.1|conserved Plasmodium protein, unknown function|No quadruplexes|0.0|
|*P. vivax*|1|EDL46376.1|hypothetical protein, conserved|Promoter|66.0|
|*P. gaboni*|1|KYO01775.1|hypothetical protein PGSY75_0607900|No quadruplexes|0.0|
|*P. knowlesi*|1|CAA9989250.1|conserved Plasmodium protein, unknown function|Promoter|69.0|
|*P. yoelii*|1|VTZ71504.1|conserved Plasmodium protein, unknown function|No quadruplexes|0.0|

![6_start](https://user-images.githubusercontent.com/60808642/174450461-634dea36-d9cf-4cd8-a6ea-741fb8c27815.png)
![6_end](https://user-images.githubusercontent.com/60808642/174450463-0a90f301-997f-4607-b0a9-c23f1eef30f8.png)
![cluster6](https://user-images.githubusercontent.com/60808642/174450095-6640e936-2826-4b38-ac33-e199467a237b.png)

#### Cluster 7
|Species|Proteins in the cluster|IDs of encoded proteins|Functions of encoded proteins|Location of quadruplexes|SCORE|
|-|-|-|-|-|-|
|*P. falciparum*|1|CZU00005.1|PhIL1 interacting protein PIP3|No quadruplexes|0.0|
|*P. vivax*|1|EDL46939.1|hypothetical protein, conserved|Promoter (outside chart)|61.0|
|*P. gaboni*|1|KYN95905.1|hypothetical protein PGSY75_1430800|No quadruplexes|0.0|
|*P. knowlesi*|1|CAA9990315.1|PhIL1 interacting protein PIP3, putative|Promoter|73.0|
|*P. yoelii*|1|VTZ78888.1|PhIL1 interacting protein PIP3, putative|No quadruplexes|0.0|

![7_start](https://user-images.githubusercontent.com/60808642/174450468-213d553f-0a0d-4fc6-bbfe-914ff2446fd6.png)
![7_end](https://user-images.githubusercontent.com/60808642/174450470-7fc74472-20af-4f28-bfa2-1d7a8c2edf57.png)
![cluster7](https://user-images.githubusercontent.com/60808642/174450103-ad311b4f-af8e-481f-84a9-1c49c0593bc8.png)

#### Cluster 8
|Species|Proteins in the cluster|IDs of encoded proteins|Functions of encoded proteins|Location of quadruplexes|SCORE|
|-|-|-|-|-|-|
|*P. falciparum*|1|CCAD51521.1|glideosome-associated protein 40, putative|No quadruplexes|0.0|
|*P. vivax*|1|EDL44073.1|hypothetical protein, conserved|Promoter|72.0|
|*P. gaboni*|1|KYO02352.1|putative glideosome-associated protein 40|No quadruplexes|0.0|
|*P. knowlesi*|1|CAA9988681.1|glideosome-associated protein 40, putative|Promoter|61.0|
|*P. yoelii*|1|VTZ79441.1|glideosome-associated protein 40, putative|No quadruplexes|0.0|

![8_start](https://user-images.githubusercontent.com/60808642/174450475-164da43b-265e-4bde-b3bd-4585221a48cb.png)
![8_end](https://user-images.githubusercontent.com/60808642/174450477-cf3d50a4-7027-49d5-9dc2-004d74d0925b.png)
![cluster8](https://user-images.githubusercontent.com/60808642/174450110-0624365a-8989-467b-88a1-1940ce5ec269.png)

#### Cluster 9
|Species|Proteins in the cluster|IDs of encoded proteins|Functions of encoded proteins|Location of quadruplexes|SCORE|
|-|-|-|-|-|-|
|*P. falciparum*|1|CAD52315.1|U4/U6.U5 tri-snRNP-associated protein 2, putative|No quadruplexes|0.0|
|*P. vivax*|1|EDL47370.1|ubiquitin carboxyl-terminal hydrolase, putative|Gene|66.0|
|*P. gaboni*|1|KYN96971.1|putative U4/U6.U5 tri-snRNP-associated protein 2|No quadruplexes|0.0|
|*P. knowlesi*|1|CAA9990777.1|U4/U6.U5 tri-snRNP-associated protein 2, putative|Gene|66.0|
|*P. yoelii*|1|VTZ81486.1|U4/U6.U5 tri-snRNP-associated protein 2, putative|No quadruplexes|0.0|

![9_start](https://user-images.githubusercontent.com/60808642/174450484-aef3618d-3e5a-40c3-b36c-083b961f3e02.png)
![9_end](https://user-images.githubusercontent.com/60808642/174450486-8e724db2-469f-4f11-b8d9-064c23714afb.png)
![cluster9](https://user-images.githubusercontent.com/60808642/174450114-b2a2f071-1268-4406-b105-a612e30771b8.png)

#### Cluster 10
|Species|Proteins in the cluster|IDs of encoded proteins|Functions of encoded proteins|Location of quadruplexes|SCORE|
|-|-|-|-|-|-|
|*P. falciparum*|1|CAX64058.1|NIMA related kinase 4|No quadruplexes|0.0|
|*P. vivax*|1|EDL43099.1|serine/threonine-protein kinase NEK4, putative|Gene|73.0|
|*P. gaboni*|1|KYO01330.1|NIMA related kinase 4|No quadruplexes|0.0|
|*P. knowlesi*|1|CAA9986571.1|NIMA related kinase 4, putative|Gene|57.0|
|*P. yoelii*|1|VTZ75464.1|NIMA related kinase 4|No quadruplexes|0.0|

![10_start](https://user-images.githubusercontent.com/60808642/174450491-6e1b7a00-943d-444e-b8d1-aa407d3ff280.png)
![10_end](https://user-images.githubusercontent.com/60808642/174450495-5c257b26-d17c-4155-b3a5-153c0293f33a.png)
![cluster10](https://user-images.githubusercontent.com/60808642/174450119-5381a020-df55-4434-afa8-2d3f3ec83c0f.png)
