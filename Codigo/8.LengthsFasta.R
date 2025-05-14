############################################################################
# Obtener longitud de los reads obtenidos de nuestro ensamble de PsAu.fasta
############################################################################
install.packages("seqinr")

library(seqinr)

seqs <- read.fasta(file = "/home/acanedo/Escritorio/PracticaR_GC2025/contigs.fa")

sequence_lengths <- sapply(seqs, function(x) {
  length(x)
})
lengthSeq <- (getLength(seqs))

hist(sequence_lengths, main = "Sequence Length Distribution", xlab = "Sequence Length", ylab = "Frequency")
