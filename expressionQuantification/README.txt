The current pipeline involves:

python readsCleaner.py				####	for both RPF and mRNA files
python genomeIndexer.py			###		run just once to create genome index for STAR
python readsMapper_mRNA.py		### 	mapping mRNA reads into genome
python readsMapper_RPF.py			###	    	mapping RPF reads into genome
python readsCounter_mRNA.py		###	    	counting reads from mRNA samples
python readsCounter_RPF.py			###	    	counting reads from RPF samples
