
import os,sys

trunk='/Volumes/omics4tb2/alomana/projects/TLR/data/ecoli_GSE53767/'
folders=os.listdir(trunk)
wfolders=[folder for folder in folders if 'counts' in folder]

print(wfolders)
wfolders.sort()

for folder in  wfolders:
    if folder != 'counts':
        print(folder)

        files=os.listdir(trunk+folder)
        wfiles=[element for element in files if '.txt' in element]
    
        print(wfiles)
        wfiles.sort()

        for wf in wfiles:
            mapped_reads=0
            with open(trunk+folder+'/'+wf,'r') as f:
                for line in f:
                    v=line.split('\t')
                    if 'transcr' in v[0]:
                        value=int(v[1])
                        mapped_reads=mapped_reads+value
            print(wf,int(mapped_reads/1e3))
    print('------')
        
