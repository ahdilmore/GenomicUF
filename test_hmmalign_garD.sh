hmmfetch Pfam-A.hmm GD_AH_C > garD_GD_AH_C.hmm
hmmalign -o garD_GD_AH_C_test.msa --trim garD_GD_AH_C.hmm garD.fa
hmmfetch Pfam-A.hmm SAF > garD_SAF.hmm
hmmalign -o garD_SAF_test.msa --trim garD_SAF.hmm garD.fa  
