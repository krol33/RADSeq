# RADSeq
Stacks wrapper

date : 29/12/15
updated : 01/03/2016
auteur : Maria Bernard
maria.bernard@jouy.inra.fr


si vous n'avez pas de SVN local il vous suffit simplement de créer un répertoire
mkdir /save/USER/git
cd /save/USER/git

git clone https://github.com/mariabernard/RADSeq.git

puis faites en une copie dans votre save au niveau d'un dossier dédié à votre projet. Les scripts sont à modifier, attention de ne pas modifier les versions du dépôt, vous n'avez pas les droits en écriture mais en faisant ça vous perdrez à chaque mise à jour du dépôt ou à chaque projet RADSeq les commandes que vous avez lancé pour un projet en particulier.
mkdir /save/USER/MON_PROJET_RAD
cd /save/USER/MON_PROJET_RAD
cp /save/USER/SVN_forge_DGA/RADSeq/* .

				

## 				1) LES DONNEES D'ENTREE

Vous avez 2 fichiers à créer

### a) indiv_barcode.txt
Ce fichier est un fichier tabulé qui décrit les individus en 4 colonnes:

		INDEX	NAME	RUN_NAME	BARCODE_SEQUENCE
		1	EFFICACE_113	EFFICACE_Run1	GCTATGTGC
		2	EFFICACE_43	EFFICACE_Run1	CTAGACGCG
		3	EFFICACE_209	EFFICACE_Run2	TGACTATAG	

- L'INDEX [OBLIGATOIRE] servira dans tous les fichiers de résultats de Stacks à identifier vos individus (généralement la deuxième colonne des fichiers).
- NAME [OBLIGATOIRE] est le nom des individus qui ont été séquencé, seul condition tous les noms doivent être unique. S'il y a des réplicats, ils doivent être identifiés différement d'une manière ou d'une autre. Biensûre pas de caractère spéciaux ou d'espacement dans les noms (ceci est valable tout le temps pour n'importe quel nom de variable ou de fichier)
	- RUN_NAME [OBLIGATOIRE] est le nom du run dans lequel l'échantillon a été séquencé. Ce nom est l'équivalent du préfix des fichiers fastq stockés dans le dossier data. Le premier script qui permet de faire le démultiplexage ira chercher automatique dans le dossier data les fichiers qui commencent par le nom du run fourni dans indiv_barcode.txt et suivi de "_1.fq.gz" et "_2.fq.gz" . Le RUN_NAME servira également à construire l'arborescence des dossiers ensuite
	- BARCODE_SEQUENCE [OPTIONNEL] est la séquence qui permet d'identifier l'individus dans le run. Si vos données sont démultiplexées alors vous n'avez pas besoins de fournir cette information et vous pouvez vous reporter à l'arborescence de fichier indiqué en 1)b)1).


### b) population.map
Ce fichier est également un fichier tabulé à 2 colonnes:

		NAME	POP_INDEX
		EFFICACE_43	1	
		EFFICACE_113	2
		EFFICACE_209	2

- NAME reprend les noms des individus indiqué dans indiv_barcode.txt
- POP_INDEX	identifie par un nombre une population ou un groupe d'individus. Si vous êtes dans le cas de données de familles, vos populations correspondent aux générations de la famille par exemple. Vous pouvez également avoir un groupe d'individus controle (généralement des individus happloïdes doublés). Si vous avez une population controle indiqué la en population "1".

Ce fichier sert à différentes étapes. Il est nativement nécessaire pour les analyses de population mais je l'utilise également pour identifier une population controle qui permet de filtrer les génotypes, et/ou pour faire des paramètrages spécifiques en fonction de groupe d'individus (haploïde doublé ou non). Ces étapes de filtres et de paramètrage doivent être adaptées à votre projet. Il permet également d'ordonner les individus dans le tableau final de génotypage. En effet, on décrit souvent indiv_barcode.txt en fonction des run, alors que population map est plutôt écrit en fonction des populations/générations, et donc c'est un ordre plus logique dans le tableau final. 

## Structure du répertoire de travail
Créez un dossier data (pour les données brutes) et un dossier SGE_out (pour les log)
	
Nommé par défaut data (l'important est que ce soit cohérent avec ce qui est indiqué dans le script), contient simplement les fichiers fastq pairés compressés du séquençage.
Le nom de fichier est parfois long, et la nomenclature de ce nom dépend des plateforme de séquençage. Je vous conseilles donc de faire des liens des fichiers brutes vers le dossier data et de renommer les liens pour avoir quelque chose de plus manipulable. 

!!! Attention toujours en respectant RUN_NAME_1.fq.gz et RUN_NAME_2.fq.gz !!!

Le nom des runs doit correspondre à ce que vous avez indiqué dans indiv_barcode.txt


## 				2) ADAPTATION DES SCRIPTS
Le pipeline est constitué de 5 scripts bash correspondant à 5 étapes clés plus deux étapes optionnelles:

		- stacks_step0.sh: preprocessing des données
		- stacks_step1.sh: clustering individuel des lectures 1
		- stacks_step2.sh: construction d'un catalogue de locus putatifs
		- stacks_step3.sh: mapping des individus sur le catalogue
		- stacks_step4.sh: genotypage/haplotypage des individus par rapport au catalogue et statistiques génétiques
		- stacks_step5.sh (optionnelle) : filtre des génotypes et tableaux de statistiques allèliques/genotypiques
		- stacks_step6.sh (optionnelle) : assemblage des mini contig

### a) généralité valables pour tous les scripts
Certaines modifications sont communes à tous les scripts.:

		- rechercher remplacer "[Project_path_dir]" par le chemin absolu qui mène à votre répertoire de travail sur votre work dédié au projet
		- rechercher remplacer [Proj_Name] par le nom de votre projet, il servira à nommer vos jobs sur le cluster et leurs fichiers de sortie
		- rechercher remplacer [Script_path_dir] par le chemin absolu du répertoire ou vous avez fait votre copie dédiée au projet du dépôt SVN qui contient tout vos scripts
		- rechercher remplacer [libR_path_dir] par le chemin absolu du répertoire du dépot libR
	
Les premières sections de scripts définissent les variables de types chemin des dossiers:
	
	- vérifier que ces chemins sont bien corrects chez vous. Normalement si vous respectez l'arborescence proposée (dossier des données = data, dossier des log = SGE_out , INDIV_FILE = indiv_bacode.txt ....) vous n'avez rien à changer.


### b) stacks_step0.sh: preprocessing des données
Une fois le script adapté avec les généralités, vous pouvez modifier les variables dédiées au programme:

		- Enzyme : le(s) nom(s) d'enzyme de restriction utilisée. Si double digestion inscrire les deux noms séparés par un espace
		- align : 0 ou 1 pour savoir si souhaitez ensuite travailler sur données alignées ou de novo.
		- paired: 0 ou 1 pour indiquer si vos données sont des données pairées ou single end 
		- déréplication : autrement dit recherche et suppression des duplicats PCR. Cela ne peut être utilisé que pour un RAD classique avec séquençage pairé. mettre dereplication=1
		- BMM : Le nombre de mismatch autorisé dans le barcode(0 est conseillé dans un premier temps, si vous ne récupéré pas suffisament de lectures alors essayé avec un mismatch)
		- TRIM indique le nombre de barcode: 1, seule la lecture 1 contient le barcode, 2 les deux lectures contiennent le barcode 
		- LEN (optionnel): la taille maximum que vous souhaitez conserver après trimming des barcodes. (ne rien mettre si barcode de même taille).
	
Vous n'avez plus qu'à lancer le script avec la commande suivante:
	qsub [Sript_path_dir]/stacks_step0.sh

Si vos données sont déjà démultipléxées, faite des liens de vos fichier fq ou fq.gz en respectant l'arborescence qui aurait résulté du démultiplexage autrement dit:

		dans le dossier preprocessing:
			- créer un dossier par run (les noms des runs sont indiqués dans la troisième colonne du fichier indiv_barcode.txt)
			- dans chaque dossier "run" faites des liens de vos fichiers fastq ou fastq.gz en les nommant si ce n'est pas le cas : [NOM_INDIVIDUS]_1.fq et [NOM_INDIVIDUS]_2.fq (ou avec l'extension gz si les fichiers sont compressés). Les [NOM_INDIVIDUS] sont indiqués dans la 2e colonne du fichier indiv_barcode.txt
			- dans le scripts stacks_step0.sh, commentez la partie "DEMULTIPLEXAGE" (ajouter des "#" en début de ligne)

Vous trouverez quelques statistiques bilan du preprocessing dans le fichier SGE_out/stacks0.out et dans le dossier stat/preprocessing.


### c) stacks_step1.sh: clustering individuel des lectures

Cette étape correspond à la clusterisation des lectures de chaque individus et à la détection des SNPs pour chaque individus.
	Les programmes ustacks ou pstacks est le coeur de cette étape.
	
	OPTIONS:
		- GENOME="" : chemin vers le génome de référence au format fasta. Celui ci doit déjà être indexé pour BWA. Si vous souhaitez travailler en de novo, laisser vide.
		- POP_HD="" : les identifiants de population des happloïdes doublés (s'il y en a), séparés par un espace. Ces identifiants sont ceux du fichier population.map
		- DD_RAD = 0 ou 1 pour indiquer si vos données sont digérées 1 ou deux fois 
		- MIN_DEPTH=3 : la couverture minimale d'un stacks autrement dit d'un allèle (-m dans ustacks) ou d'un locus (-m dans pstacks)
		- MAX_PRIM_DIST=2 : le nombre de mismatch autorisés entre stacks, autrement dit le nombre potentiel de SNP sur un locus (-M dans ustacks)
		- MAX_SEC_DIST=4 : le nombre de mismatch autorisés pour agréger des "secondary reads", augmenter la couverture des locus (-N par défaut = -M+2 dans stacks)
		- STACKS_OPT="" (optionnel): si vous souhaitez ajouter d'autres option ustacks (sachant qu'il y a déjà -r et -d ). Gardez les "" .
		
Différence entre happloïde doublé ou non:
		- le nombre maximal de stacks par cluster, autrement dit le nombre d'allèles possibles par locus (--max_locus_stacks) : pour les individus happloïdes doublés tous les locus sont attendus sans SNP donc avec un allèle + une marge d'erreur --max_locus_stacks = 2, suivant la même théorie pour un individus hétérozygote 2 allèles + erreur = --max_locus_stacks = 3 . Cette option est indiquée en dure dans les commande ustacks
			
Les individus happloïdes doublés permettent de contrôler la clusterisation. Dans tous les cas on attend plus de site homozygote qu'hétérozygote, et dans le cas des happloïdes doublés la différence de proportions doivent être encore plus marquée.

Les paramètres indiqués sont ceux que j'ai utilisé jusque là. A vous de les adapter en fonction de la qualité de vos données (notamment la couverture des allèles : -m ) et l'ajout ou non des secondary reads. 
		
Vous n'avez plus qu'à lancer le script avec la commande suivante:
	qsub [Sript_path_dir]/stacks_step1.sh	

Vous trouverez un certain nombre de statistiques descriptives (summary_statistics.txt et des graph) dans le dossier stat/stacks1

### d) stacks_step2.sh: construction d'un catalogue de locus putatifs

Cette seconde étape de l'analyse des RAD, correspond à la constitution d'un catalogue de locus à partir des clusters individuels précédents. Il s'agit de sélectionner tout ou partie des individus permettant de représenter la majeure partie des polymorphismes de l'ensemble des individus. Le programme sousjacent est le programme cstacks de la suite Stacks (http://creskolab.uoregon.edu/stacks/comp/cstacks.php )
	
Dans le cas d'analyse de famille, on sélectionne théoriquement la F0 et F1 (ou uniquement la F1). Si ces individus sont mal séquencés ou présentent des résultats de clusterisation moyens, vous pouvez également raisonner en analyse de population et prendre un certain nombre de F2 qui eux auront de bonnes couvertures et à priori une bonne clusterisation. 
	
Dans le cas d'analyse de population on sélectionne théoriquement tous les individus. Si vous avez une population contrôle type happloïdes doublés, ne l'intégrez pas au catalogue.
	
	OPTIONS: 
		- MISMATCH=1 : le nombre de mismatch autorisé entre 2 séquences consensus de cluster individuel (en plus des SNP déjà détecté dans chaque individus).Si vos populations de départ sont supposés plus éloignées (souches ou espèces différentes ) vous pouvez augmenter (un peu !) ce paramètre.
		- KNOWN_CAT= : si vous avez déjà un catalogue stacks fait avec le même protocole RAD, indiquez le chemin complet vers les résultats de sortie. (ex : /PATH/VERS/CSTACKS_DIRECTORY/batch_1.catalog )
		- align = 0 ou 1 pour indiquer si vous avez travaillé sur données alignées (pstacks) ou sur un clustering de novo (ustacks)

		# INDIVIDUALS SELECTION (1 option au choix, si aucune option tous les individus sont pris en compte). Indiquez séparé par des espaces les index de pop ou les individus que vous souhaitez agrégés dans le catalogue ou au contraire ne pas mettre. Si vous souhaitez tout prendre laissez vide.
			POP_IN=""
			POP_OUT=""
			INDIV_IN=""
			INDIV_OUT=""

		# VARIABLE CLUSTER option de cpu et mémoire. Normalement une telle configuration doit convenir pour tous.
		MEM="16G"
		VMEM="24G"
		THREADS=10
		Attention si vous souhaitez augmenter ces paramètres, pensez bien que la réservation mémoire se fait par CPU donc ici 10*16 = 160G et 10*24 = 240G

Un fois adapté à vos données, lancez : 
	qsub [Sript_path_dir]/stacks_step2.sh	

Vous trouverez un certain nombre de statistiques descriptives (summary_statistics.txt et des graph) dans le dossier stat/cstacks


### e) stacks_step3.sh: mapping des individus sur le catalogue

Après avoir fait une clusterisation par individus et l'intersection de ces clusters dans un catalogue, il s'agit de refaire correspondre chaque cluster individuel à un locus du catalogue grâce au programme sstacks de la suite Stacks (http://creskolab.uoregon.edu/stacks/comp/sstacks.php ).
	
	OPTION:
		- seuil_locus=0 : qui n'est qu'un seuil maximum pour le calcul du nombre d'individus ayant une correspondance pour moins de S locus. Ce nombre sera également calculé pour S=10, 100, 1000, et 10000. Si vous n'avez aucune idée du nombre de locus attendus laissez ce seuil à 0
		- align = 0 ou 1 pour indiquer si vous avez travaillé sur données alignées (pstacks) ou non (ustacks)
	
Puis lancez :
		qsub [Sript_path_dir]/stacks_step3.sh	
	
NB: les individus n'ayant aucune correspondance avec le catalogue seront indiqués dans le repertoire de sortie dans le fichier indiv_to_be_removed.txt et seront potentiellement exclu de la suite de l'analyse.
	
Vous trouverez un certain nombre de statistiques descriptives (summary_statistics.txt et fichier de comptage) dans le dossier stat/sstacks

### f) stacks_step4.sh: genotypage/haplotypage des individus par rapport au catalogue et statistiques génétiques

Cette étape est la dernière du pipeline Stacks. Selon le types d'analyse elle fait appel au programme populations (http://creskolab.uoregon.edu/stacks/comp/populations.php ) ou au programme genotypes (http://creskolab.uoregon.edu/stacks/comp/genotypes.php ). 
Le scripts stacks_step4.sh, va lancer ce dernier programme puis réordonner le tableau d'haplotypage en fonction de l'ordre indiqué dans le fichier population.map, et filtre les locus monomorphe des locus polymorphes (qui vous intéressent)
	
	OPTIONS:
		- ANALYSIS_TYPE="populations" : "populations" ou "genotypes" (pour l'analyse de famille)
		- seuil_locus=0 : qui correspond au nombre minimum de locus génotypés par individus que vous souhaitez pour identifier des individus de mauvaise qualité et donc les filtrer dans la prochaine étape (mettez 0 ou rien si vous ne souhaitez pas filtrer d'individus)
		
Puis lancez:
	qsub [Sript_path_dir]/stacks_step4.sh	

Vous trouverez un certain nombre de statistiques descriptives (summary_statistics.txt et fichier de comptage) dans le dossier stat/populations ou stat/genotypes


### g) stacks_step5.sh (optionnelle) : filtre des génotypes et tableaux de statistiques allèliques/genotypiques

Ce script permet de filtrer les locus en fonction d'une population controle type happloïde doublée. 
	
	OPTIONS:
		- POP_CONTROLE="". identifiants des populations haploïdes doublés
		- NB_INDIV=1. nombre minimum de de genotype Hz dans les pop control pour déclarer le locus comme mauvais
		
Lancez ensuite la commande 
	qsub [Sript_path_dir]/stacks_step5.sh	


### h) stacks_step6.sh (optionnelle) : assemblage des mini contig	

Cette dernière étape permet d'assembler les lectures d'une sélection de locus d'intérêt en contig puis de refaire une détection de SNP classique via GATK pour valider les SNP détectés avec Stacks et de trouver la localisation potentielle de ces locus sur le génome de référence.
	
Vous devez donc adapter ce script en fonction de vos données.
	
	OPTIONS:
		- ENZ=""
		- QUAL=		# filtre qualitatif des SNP détecté par GATK: nb indiv * 50 (Qualité moyenne > 50 (sur une echelle de 0 à 100)
		- AN=		# filtre quantitatif des SNP détecté par GATK: 50% des individus génotypés = nombre d'individus
		- GENOME= 	# fichier fasta du génome de référence si possible

		# INPUT
		- SAMPLE_READS_DIR=	# dossier contenant les lectures R1 et R2 apres process_radtag de chaque individus (dossier checkRAD)
		- READS_TYPE="gzfastq"	# fastq ou gzfastq
		- STACKS_RES_DIR=		# dossier de résultat stacks (normalement populations/input ou genotypes/input)
		- WHITE_LOCUS_LIST=		# fichier des locus à assembler (un identifiant de locus par ligne)
	
Ce programme peut être très long. Je vous conseille de lancer le script avec un fichier contenant 2 locus d'interet pour vérifier que tout se passe bien avant de le relancer sur votre liste complète de locus d'intérêt.
		
lancez qsub [Sript_path_dir]/stacks_step6.sh	
	
	

