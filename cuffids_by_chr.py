Type "help", "copyright", "credits" or "license" for more information.
>>> lines=open('/Volumes/KemiTegel/ProstateTumorAnalysis/DataAnalysisR/cuffannotation_rRNAd/merged_requireStrand.gtf').readlines()
>>> len(lines)
3219489
>>> lines2=[line.split('\t') for line in lines]
>>> lines2=[line.split('\t') for line in lines if line.split('\t')[0] == '2']
>>> len(lines2)
245152
>>> trans_ids = [line[8].split('; ')[1][15:-1] for line in lines2]
>>> trans_ids[0]
'TCONS_00210065'
>>> trans_ids2 = set(trans_ids)
>>> len(trans_ids)
245152
>>> len(trans_ids2)
29869
>>> out=open('/Volumes/KemiTegel/ProstateTumorAnalysis/DataAnalysisR/cuffannotation_rRNAd/chr2_transcript_ids.txt','w')
>>> for l in trans_ids2:
...     out.write(l)
...
>>> out.close()
>>> out=open('/Volumes/KemiTegel/ProstateTumorAnalysis/DataAnalysisR/cuffannotation_rRNAd/chr2_transcript_ids.txt','w')
>>> for l in trans_ids2:
...     out.write(l + '\n')
...
>>> trans_ids2[-1]
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
TypeError: 'set' object does not support indexing
>>> list(trans_ids2)[-1]
'TCONS_00216992'
>>> for l in trans_ids2:
...
KeyboardInterrupt
>>> out.close()
