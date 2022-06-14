## Поиск участков Z-ДНК в промоторах генов представителей рода Plasmodium (тип Apicomplexans)
Для анализа генома на предмет наличия Z-ДНК был выбран тип *Apicomplexa*, род *Plasmodium*.

Данные по видам *P. falciparum*, *P. vivax*, *P. knowlesi*, *P. gaboni*, *P. yoelii* (INSDC) были скачаны с https://www.ncbi.nlm.nih.gov/genome/browse#!/eukaryotes/refseq_category:representative и https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA. Причём файлы .fasta и .gb были скачаны из браузера NCBI с помощью скрипта на bash, по одному файлу на хромосому: итого 70 .fasta файлов и 70 .gb файлов. Файлы же feature_table.txt были скачаны через ftp.

|Вид|Число хромосом|Длина генома|Количество анн. генов|Доля анн. генов на геном|Число участков со значимым ZH_SCORE*|Длина участков со значимым ZH-SCORE*|
|-|-|-|-|-|-|-|
|*P. falciparum*|14**||||||
|*P. vivax*|14||||||
|*P. gaboni*|14||||||
|*P. knowlesi*|14**||||||
|*P. yoelii*|14**||||||

*А именно с ZH-SCORE >= 500.

** Помимо 14 хромосом имелись также данные по митохондриальному геному и геному апикопласта, однако там нашёлся всего один участок с ZH-SCORE >= 500, поэтому их не рассматривали.

#### Кластер 1
|Вид|Генов в кластере|ID кодируемых белков|Функция кодируемых белков|Расположение Z-ДНК|ZH_SCORE|
|-|-|-|-|-|-|
|*P. falciparum*|1|CZT99390.1|cAMP-dependent protein kinase regulatory subunit|||
|*P. vivax*|1|EDL47592.1|cAMP-dependent protein kinase regulatory subunit, putative|||
|*P. gaboni*|1|KYN98034.1|cAMP-dependent protein kinase regulatory subunit|||
|*P. knowlesi*|1|CAA9991018.1|cAMP-dependent protein kinase regulatory subunit, putative|||
|*P. yoelii*|1|VTZ81713.1|cAMP-dependent protein kinase regulatory subunit, putative|||


#### Кластер 2
|Вид|Генов в кластере|Гены|Функция генов|Расположение Z-ДНК|ZH_SCORE|
|-|-|-|-|-|-|
|*P. falciparum*|1|CAD51193.2|zinc finger protein, putative|||
|*P. vivax*|1|EDL45264.1|D13 protein, putative|||
|*P. gaboni*|1|KYO00773.1|putative zinc finger protein|||
|*P. knowlesi*|1|CAA9986906.1|zinc finger protein, putative|||
|*P. yoelii*|1|VTZ76436.1|zinc finger protein, putative|||


#### Кластер 3
|Вид|Генов в кластере|Гены|Функция генов|Расположение Z-ДНК|ZH_SCORE|
|-|-|-|-|-|-|
|*P. falciparum*|1|CZT99041.1|conserved Plasmodium protein, unknown function|||
|*P. vivax*|1|EDL45764.1|hypothetical protein, conserved|||
|*P. gaboni*|1|KYN99013.1|hypothetical protein PGSY75_1137600|||
|*P. knowlesi*|1|CAA9988401.1|conserved Plasmodium protein, unknown function|||
|*P. yoelii*|1|VTZ78317.1|conserved Plasmodium protein, unknown function|||


#### Кластер 4
|Вид|Генов в кластере|Гены|Функция генов|Расположение Z-ДНК|ZH_SCORE|
|-|-|-|-|-|-|
|*P. falciparum*||||||
|*P. vivax*||||||
|*P. gaboni*||||||
|*P. knowlesi*||||||
|*P. yoelii*||||||

#### Кластер 5
|Вид|Генов в кластере|Гены|Функция генов|Расположение Z-ДНК|ZH_SCORE|
|-|-|-|-|-|-|
|*P. falciparum*||||||
|*P. vivax*||||||
|*P. gaboni*||||||
|*P. knowlesi*||||||
|*P. yoelii*||||||
