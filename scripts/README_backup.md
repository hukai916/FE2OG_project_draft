Pipeline: Step1-5

Step1: get a list of seeds for FE2OG structures.
  proc1 - SCOPe search
  step: search for iron and ketoglutarate dependent protein in SCOPe database.
  results: found superfamily b.82.2(Clavaminate synthase-like) that covers 38 representative domains; further filtered (remove domains with extreme lengths) to 30 domains (/Users/Kai/Desktop/GitHub/FE2OG_project_draft/SCOPe_domain).

  proc2 - Fatcat search
  step: search each of the 30 domains against 90% non-redundant PDB database using Fatcat server(http://fatcat.sanfordburnham.org/).
  results: 30 lists of hits (/Users/Kai/Desktop/GitHub/FE2OG_project_draft/Fatcat)

  proc3 - Combine and filter Fatcat results
  step: combine the 30 lists and get rid of redundant hits.
  results: 1423 unique pdb hits

  proc4 - Further filter Fatcat results
  step: check the frequencies and qualities of those 1423 unique hits.
  results: hits with less than 5 occurrences are basically fragments and false positives, therefore, keep 297 domains that have at least 5 occurrences. (/Users/Kai/Desktop/GitHub/FE2OG_project_draft/Fatcat/297_domain_5.txt)
    Alternative proc4 - use P-value and alignment-length percentage (termed Q-value) as criteria
                step: to determine the cutoff, investigate the hits for 1ds1, 2fct and 3o2g.
                For 1ds1, the minimal Q-value is 0.43 and is still FE2OG enzyme;
                for 2fct, the minimal Q-value is 0.06 (2k9iA), which is not FE2OG; extend to 0.55 (4lt5A) and found the confirmed FE2OG;
                for 3o2g, the minimal Q-value is 0.19, which is not FE2OG, extend to 0.55 (3r1j) and found the first confirmed FE2OG.
                We first try to set the Q-value cutoff to 0.50 and ended up with 329 unique hits this way. We know 329 covers all of the 297 hits and is too much. Then, try cutoff at 0.80, we get 117 hits, compare these 117 to 80 manually confirmed domains, only 6/80 is not covered by these 117. If using 0.75, 153 hits cover all of the 80 confirmed domains. We can use 0.75 at the moment.
                (ls ../Fatcat/\*list | xargs -i python parseFatcat.py {} | sort | uniq > 329.id). Therefore, try to adjust the hit cutoff: 0.60 lead to 280 hits, 0.7 lead to 175 hits and 0.8 lead to 117 hits, 0.85 lead to 88 hit.


  proc5 - Check these 297 hits and summarize. (when cutoff set to 0.5)
  results: after manual investigation, 80 are characterized FE2OG enzymes, 17 are potential FE2OG enzymes, the rest may not be FE2OG enzymes at all.
  (/Users/Kai/Desktop/GitHub/FE2OG_project_draft/Fatcat/297_domain_summary.xls)

  proc5 - choose the set of seeds to test

  Checkpoint: how many of the 19 clusters that we have analyzed before are covered by the 80 hits we found?
  -- all but H6H and LolO, our method should be valid and would identify most FE2OG enzymes in the database.

Step2: find close homologs for each of the seed from Step1
(153 using cutoff at 0.75, this 153 set is supposed to be 'real' Fe2OG enzymes that would have a structure associated with them).

proc1 - blast parameter:
        blastp maximum hit number set to 5000
proc2 - filter parameter: use cdhit to speed up the clustering process.
        install cdhit: https://github.com/weizhongli/cdhit/wiki/2.-Installation
        (xcode-select --install: use this to install command line tools which allows you to use most compilers without having to specify the location of the headers manually)
    subproc2 - add pdb sequence on the top of each blastp hits to make sure the cd-hit output will put pdb sequence on top of the cluster file as well. (ls ../Fasta/\*.fasta | xargs -n 1 -P 10 python clusterReorder.py)

Step3: run pairwise tree building and SSN. Instead, try cd-hit at different cutoffs
(Because some of the 153 sequences might be very similar and should be grouped into the same group, after blast, recluster the unique hits at certain level to see how many groups left.)
(It would be hard to automate the selection of clusters from either tree or SSN.
  Instead, we can use cd-hit to re-cluster the hits to low identity clusters. The min cutoff of cd-hit is 0.4 which gives 861 clusters for all 72526 sequences. This is still too much, we can use psi-cd-hit to cluster down to 30% or even lower to see how it would look like)(The installation of psi-cd-hit is as below, rather, try cd-hit webserver)
http://weizhongli-lab.org/cdhit_suite/cgi-bin/result.cgi?JOBID=1530311159

proc1 -- determine the correct clusters.
(How to deal with the overlapping among clusters?)
To find out the correct cluster:
Step1 -- Remove hits with duplicated protein ID: under /Users/Kai/Desktop/GitHub/FE2OG_project_draft/Blastp_cluster folder
cat all_unique_id.txt | xargs -n 1 -P32 efetch -db protein -format fasta -id > all_unique.fasta
（this would download the full length of the protein that is not preferred, therefore, use another customized script to filter out the needed unique seqs）
python filter_uniq.py ../Blastp_cluster/all_seq.fasta > ../Blastp_cluster/all_filter.fasta

proc2 -- Try clustering of the result at different cutoffs
psi-cd-hit.pl -i all_unique.fasta -o run30 -c 0.3 -para 8 -blp 4
psi-cd-hit.pl -i all_unique.fasta -o run30 -c 0.3 -para 8 -blp 4

psi-cd-hit.pl -i all_filter.fasta -o filter30 -c 0.3 -para 8 -blp 4
psi-cd-hit.pl -i all_filter.fasta -o filter25 -c 0.25 -para 8 -blp 4
psi-cd-hit.pl -i all_filter.fasta -o filter25 -c 0.2 -para 8 -blp 4

0.2 might be a proper cutoff for clustering (108 clusters, saved in filter20 folder).

proc3 -- reorder the sequences by putting the structure seq on top of each cluster.
python parse_clstr.py filter20.clstr 153.txt
python get_clstr.py filter20/108/* all_filter.


After we have the clusters: use the following to get the information of the clusters:
for x in Cluster_*; do head -n 1 $x; done > 114_id.txt (/Users/Kai/Desktop/GitHub/FE2OG_project_draft/Blastp_cluster/filter_0.95_0.35/Cluster)
for x in Cluster_*; do python findSeed.py $x; done > findSeed_result.txt
Then, combine the findSeed_result.txt into Cluster_114_info.xlsx.



Step4: find the equivalent sites by using TMalign for each clusters



Step5: Try different methods on the clusters, say, PCA, Onehot etc.
proc1 - try on entire space of clusters
proc2 - try on closely related clusters
(if key sites cluster them together while full sequence similarity tree separate them apart.)


Code:
Step1-proc4:
ls ../Fatcat/\*.list | xargs -i python parseFatcat.py {} | sort | uniq > ../Fatcat/153_0.75_domain.txt
python parseXML.py ../Fatcat/153_0.75_domain.txt

Step2-proc1
under script folder:


Step2-proc2
under script folder:
ls ../Blastp/* | xargs -n 1 -P 10 python parseBlastp.py # run about hours, a quick alternative would be to retrieve sequences from XML itself instead of retrieving again using efetch.
For fast version:
ls ../Blastp/* | xargs -n 1 -P 10 python parseBlastp_fast_version.py # can be done in minutes.
ls ../Fasta/\*.fasta | xargs -n 1 -P 10 python clusterReorder.py

/Users/Kai/Desktop/GitHub/FE2OG_project_draft/Blastp_cluster
cat \*.reorder* > all_seq.fasta




Installation of psi-cd-hit:
Already integrated inside cd-hit package. Only need to add to PATH.
