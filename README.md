## Поиск участков Z-ДНК в промоторах генов представителей рода Plasmodium (тип Apicomplexans)
Для анализа генома на предмет наличия Z-ДНК был выбран тип *Apicomplexa*, род *Plasmodium*.

Данные по видам *P. falciparum*, *P. vivax*, *P. knowlesi*, *P. gaboni*, *P. yoelii* (INSDC) были скачаны с https://www.ncbi.nlm.nih.gov/genome/browse#!/eukaryotes/refseq_category:representative и https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA. Причём файлы .fasta и .gb были скачаны из браузера NCBI с помощью скрипта на bash, по одному файлу на хромосому: итого 70 .fasta файлов и 70 .gb файлов. Файлы же feature_table.txt были скачаны через ftp.

|Вид|Число хромосом|Длина генома|GC состав|Количество анн. генов|Доля анн. генов на геном|Число участков со значимым ZH_SCORE*|Длина участков со значимым ZH-SCORE*|
|-|-|-|-|-|-|-|-|
|*P. falciparum*|14**|23 292 622|19.34|5387|0.59|1064|13 714|
|*P. vivax*|14|22 621 071|44.87|5389|0.62|37018|476 891|
|*P. gaboni*|14|18 878 758|17.84|5198|0.65|617|8 162|
|*P. knowlesi*|14**|23 938 832|38.61|5326|0.59|18210|234 759|
|*P. yoelii*|14**|23 043 114|21.74|6039|0.56|1844|24 712|

*А именно с ZH-SCORE >= 500.

** Помимо 14 хромосом имелись также данные по митохондриальному геному и геному апикопласта, однако там нашёлся всего один участок с ZH-SCORE >= 500, поэтому их не рассматривали.

#### Для 5 видов proteinortho построил 5319 кластеров.

Число кластеров, состоящих из n геномов имеет средующее распределение:

![clust_by_number_of_genoms](https://user-images.githubusercontent.com/60808642/173706235-3e63ee97-388a-4c93-91a8-9709d930ce60.png)


Число кластеров, в которые входят n белковых последовательностей, имеет слдующее распределение:

![clust_by_number_of_proteins](https://user-images.githubusercontent.com/60808642/173706248-663fa365-f091-4f64-818e-d7a863e28d49.png)


#### Кластер 1
|Вид|Генов в кластере|ID кодируемых белков|Функция кодируемых белков|Расположение Z-ДНК|ZH_SCORE|
|-|-|-|-|-|-|
|*P. falciparum*|1|CZT99390.1|cAMP-dependent protein kinase regulatory subunit||583.4285|
|*P. vivax*|1|EDL47592.1|cAMP-dependent protein kinase regulatory subunit, putative||138924.1|
|*P. gaboni*|1|KYN98034.1|cAMP-dependent protein kinase regulatory subunit||583.4285|
|*P. knowlesi*|1|CAA9991018.1|cAMP-dependent protein kinase regulatory subunit, putative||583.4285|
|*P. yoelii*|1|VTZ81713.1|cAMP-dependent protein kinase regulatory subunit, putative||583.4285|

![cluster1](https://user-images.githubusercontent.com/60808642/174160532-0c6192f2-66b1-4b3b-a79f-9f6f77b8811c.png)

#### Кластер 2
|Вид|Генов в кластере|ID кодируемых белков|Функция кодируемых белков|Расположение Z-ДНК|ZH_SCORE|
|-|-|-|-|-|-|
|*P. falciparum*|1|CAD51935.1|protein phosphatase-beta||0.0|
|*P. vivax*|1|EDL45060.1|protein phosphatase-beta, putative||94590.41|
|*P. gaboni*|1|KYO00236.1|protein phosphatase-beta||705.4245|
|*P. knowlesi*|1|CAA9987523.1|protein phosphatase-beta, putative||650.9198|
|*P. yoelii*|1|VTZ77586.1|protein phosphatase-beta, putative||583.4285|

![cluster2](https://user-images.githubusercontent.com/60808642/174160560-c4cf8210-00d5-47ae-ae28-40608a534ed4.png)

#### Кластер 3
|Вид|Генов в кластере|ID кодируемых белков|Функция кодируемых белков|Расположение Z-ДНК|ZH_SCORE|
|-|-|-|-|-|-|
|*P. falciparum*|1|CAB11145.2|formate-nitrite transporter||2091.083|
|*P. vivax*|1|EDL44820.1|transporter, putative||38833.58|
|*P. gaboni*|1|KYO03279.1|putative formate-nitrite transporter||6565.992|
|*P. knowlesi*|1|CAA9987895.1|formate-nitrite transporter, putative||28780.5|
|*P. yoelii*|1|VTZ73436.1|formate-nitrite transporter, putative||0.0|

![cluster3](https://user-images.githubusercontent.com/60808642/174160577-2a079fbd-54cf-4932-a4d4-49fbbf84ca5a.png)

#### Кластер 4
|Вид|Генов в кластере|ID кодируемых белков|Функция кодируемых белков|Расположение Z-ДНК|ZH_SCORE|
|-|-|-|-|-|-|
|*P. falciparum*|1|CZT98040.1|octaprenyl pyrophosphate synthase||1743.107|
|*P. vivax*|1|EDL43227.1|prenyl transferase, putative||16770.32|
|*P. gaboni*|1|KYO03461.1|octaprenyl pyrophosphate synthase||0.0|
|*P. knowlesi*|1|CAA9986822.1|octaprenyl pyrophosphate synthase, putative||2180.071|
|*P. yoelii*|1|VTZ72360.1|octaprenyl pyrophosphate synthase, putative||583.4285|

![cluster4](https://user-images.githubusercontent.com/60808642/174160605-0062db79-4ad6-49ff-9265-e3f7d19497b8.png)

#### Кластер 5
|Вид|Генов в кластере|ID кодируемых белков|Функция кодируемых белков|Расположение Z-ДНК|ZH_SCORE|
|-|-|-|-|-|-|
|*P. falciparum*|1|CAD52619.1|peptidase, putative||612.3848|
|*P. vivax*|1|EDL44548.1|hypothetical protein, conserved||8323.257|
|*P. gaboni*|1|KYN97268.1|putative membrane protein||612.3848|
|*P. knowlesi*|1|CAA9989837.1|peptidase, putative||8323.257|
|*P. yoelii*|1|VTZ81050.1|peptidase, putative||0.0|

![cluster5](https://user-images.githubusercontent.com/60808642/174160640-7beb2aa8-bd5f-4387-9bb1-960a76760d33.png)

#### Кластер 6
|Вид|Генов в кластере|ID кодируемых белков|Функция кодируемых белков|Расположение Z-ДНК|ZH_SCORE|
|-|-|-|-|-|-|
|*P. falciparum*|1|CZT98598.1|conserved protein, unknown function||915.9191|
|*P. vivax*|1|EDL42405.1|hypothetical protein PVX_111020||13713.99|
|*P. gaboni*|1|KYN99721.1|hypothetical protein PGSY75_1034100||915.9191|
|*P. knowlesi*|1|CAA9987243.1|conserved protein, unknown function||908.3955|
|*P. yoelii*|1|VTZ74384.1|conserved protein, unknown function||0.0|

![cluster6](https://user-images.githubusercontent.com/60808642/174160667-bd6e244a-7c81-4557-8677-39a5a4bdf0b4.png)

#### Кластер 7
|Вид|Генов в кластере|ID кодируемых белков|Функция кодируемых белков|Расположение Z-ДНК|ZH_SCORE|
|-|-|-|-|-|-|
|*P. falciparum*|1|CZT98374.1|hypoxanthine-guanine phosphoribosyltransferase||583.4285|
|*P. vivax*|1|EDL44708.1|hypoxanthine phosphoribosyltransferase, putative||2276.8|
|*P. gaboni*|1|KYN99499.1|hypoxanthine-guanine phosphoribosyltransferase||0.0|
|*P. knowlesi*|1|CAA9987765.1|hypoxanthine-guanine phosphoribosyltransferase, putative||6408.923|
|*P. yoelii*|1|VTZ79941.1|hypoxanthine-guanine phosphoribosyltransferase, putative||2091.083|

![cluster7](https://user-images.githubusercontent.com/60808642/174160691-a5199f69-0862-40c3-9eb1-10568d1bc7f5.png)

#### Кластер 8
|Вид|Генов в кластере|ID кодируемых белков|Функция кодируемых белков|Расположение Z-ДНК|ZH_SCORE|
|-|-|-|-|-|-|
|*P. falciparum*|1|CZT98084.1|conserved Plasmodium protein, unknown function||2183.574|
|*P. vivax*|1|EDL43268.1|hypothetical protein, conserved||1820.585|
|*P. gaboni*|1|KYO03505.1|hypothetical protein PGSY75_0207100||2183.574|
|*P. knowlesi*|1|CAA9986780.1|conserved Plasmodium protein, unknown function||4831.423|
|*P. yoelii*|1|VTZ72434.1|conserved Plasmodium protein, unknown function||0.0|

![cluster8](https://user-images.githubusercontent.com/60808642/174160707-d2d6a252-a5bf-48c4-9459-0d02bda79429.png)

#### Кластер 9
|Вид|Генов в кластере|ID кодируемых белков|Функция кодируемых белков|Расположение Z-ДНК|ZH_SCORE|
|-|-|-|-|-|-|
|*P. falciparum*|1|CZT98918.1|conserved protein, unknown function||731.2843|
|*P. vivax*|1|EDL45645.1|hypothetical protein PVX_091900||3661.11|
|*P. gaboni*|1|KYN93456.1|hypothetical protein PGSY75_0013800||0.0|
|*P. knowlesi*|1|CAA9988282.1|conserved protein, unknown function||3661.11|
|*P. yoelii*|1|VTZ78436.1|conserved protein, unknown function||2173.083|

![cluster9](https://user-images.githubusercontent.com/60808642/174160723-50ba630e-251d-42d7-b499-bf1aba398401.png)

#### Кластер 10
|Вид|Генов в кластере|ID кодируемых белков|Функция кодируемых белков|Расположение Z-ДНК|ZH_SCORE|
|-|-|-|-|-|-|
|*P. falciparum*|1|CAD51193.2|zinc finger protein, putative||2183.574|
|*P. vivax*|1|EDL45264.1|D13 protein, putative||826.8355|
|*P. gaboni*|1|KYO00773.1|putative zinc finger protein||2183.574|
|*P. knowlesi*|1|CAA9986906.1|zinc finger protein, putative||826.8355|
|*P. yoelii*|1|VTZ76436.1|zinc finger protein, putative||2183.574|

![cluster10](https://user-images.githubusercontent.com/60808642/174160733-f5616104-11ce-4e09-a292-f43eee7be170.png)
