# HybridSystemModels
This repository contains software for Hybrid System based gene regulatory network Models (HSM). 
It consists of four tools that can be applied in sequence. The input is ... 

We provide two examples, the first is HSM lambda phage model located in file ...
![](./assets/LPH2.png)

The second is HK022 phage model located in file...
![](./assets/HK022.png)

The ’core’ models are shown in red and contain genes (cI, cII,cro) and binding sites (bOR, bcII-1) that are involved in regulatory feedback.
## HSMFeatureExtractor 

## HSMFrameRefinementExplorer

To run it use the following command line: 
```sh
python3 print_statistics.py data_folder
```
The parameter data_folder is the folder name containing state space files output by HSMSpaceAnalyzer
 
## HSMModelConverter

## HSMSpaceAnalyzer
Saliku sajā "savu daļu". Izejas faili iekš src, kompilētie gan bin, gan līmeni zemāk, no kura var palaist tā lai atpazīst lokālos pathus uz failiem (plānoju pielikty vēl īsu readme).
Folderī HSM_Models - modeļu ieejas faili, t.i. tie, kurus vajadzētu iegūt no Dārtas programmas (atskaitot hsm_loops_2.txt, par kuru pieņemam,. ka rakstīts tikai "tehniskajā
formātā"). Ģenerētie modeļu ss faili nonāk HSM_StateSpaces. Paliekam pie tiewm 9 (8 bez loops) modeļiem, kas šeit ielikti. Pāris (dažas ne "pārak smukas" atrunas):

Modeļu faili tiek padoti kā parametrs, var saukties kā pagadās, bet svarīgs ir pirmais teksta strings (ModelName) šajos failos e.g. "Lambda_Complete", tajā skaitā:

1) Programma šobrīd darbojas tikai zināmiem ModelName- resp., tiem 9, kas šeit parādās
2) Uzģenerētie telpu faili tiks saglabāt HSM_StateSpaces/[ModelName] kuram tur ir jābūt
2) Dārtai ģenerējot modeļus tie acīmredzot uzrodas no "readable" ieejas failu nosaukumiem - e.g.  Lambda_Complete.txt dod "Lambda_Complete" model name. Varam pie šī palikt, būtu gan labāk, ja ārtas prtogrammai šo ModelName varētu padot kā parametru.


