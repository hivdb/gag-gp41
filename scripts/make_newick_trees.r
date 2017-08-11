#! /usr/bin/env Rscript
library(ape)

BASEDIR = normalizePath('..')

makeTree = function(gene, category) {
  fasta = read.FASTA(sprintf('%s/data/fasta/%s%s.aln.fasta.txt', BASEDIR, gene, category))
  distResult = dist.dna(fasta, model = 'raw', pairwise.deletion = TRUE, as.matrix = TRUE)
  njResult = nj(distResult)
  write.tree(njResult, sprintf('%s/result_data/%s%s.tre', BASEDIR, gene, category))
  cat(sprintf('%s/result_data/%s%s.tre created\n', BASEDIR, gene, category))
}

makeTree('Gag', 'NNRTIs')
makeTree('Gp41', 'NNRTIs')
makeTree('PR', 'NNRTIs')
makeTree('Gag', 'PIs')
makeTree('Gp41', 'PIs')
makeTree('PR', 'PIs')