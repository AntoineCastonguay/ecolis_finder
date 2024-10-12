ecoli <- read.table("~/Documents/output3/2_result/output.txt", header = T)
ecoli_BW25113 <- read.table("~/Documents/Ecoli_R/All-genes-of-E.-coli-K-12-substr.-BW25113.txt", header = T)

colnames(ecoli) <- c('gene', 'flag', 'first_pos', 'second_pos', 'length', 'quality', 'essentiel')
# second_pos = position de la pairs

ecoli$first_pos <- ifelse(ecoli$flag == 99, ecoli$first_pos + 47,
       ifelse(ecoli$flag == 147, ecoli$first_pos + 20, ecoli$first_pos))

ecoli$second_pos <- ifelse(ecoli$flag == 99, ecoli$second_pos + 20,
                          ifelse(ecoli$flag == 147, ecoli$second_pos + 47, ecoli$second_pos))

huge_gene <- subset(ecoli,abs(ecoli$length) > 100000)
# Il y a 10 gene dont la longueur est aberente dont 3 presume essentiel.

other_flag <- subset(ecoli, ecoli$flag %in% c(65, 129, 81, 161, 97, 145, 113, 177))
good_flag <- subset(ecoli, ecoli$flag %in% c(83,99,147,163))
num_flag <- table(ecoli$flag) 
num_flag

# Les flag
# Mapped within the insert size and in correct orientation
# 99 = read paired, read mapped in proper pair, *mate* reverse strand et *first* in pair.
# 147 = read paired, read mapped in proper pair, *read* reverse strand et *second* in pair.
# 83 = read paired, read mapped in proper pair, *read* reverse strand et *first* in pair.
# 163 = read paired, read mapped in proper pair, *mate* reverse strand et *second* in pair.
# Mapped uniquely, but with wrong insert size
# 65 = read paired et first in pair
# 129 = read paired et second in pair
# 81 = read paired, read reverse strand et first in pair
# 161 = read paired, mate reverse strand et second in pair.
# 97 = read paired, mate reverse strand et first in pair.
# 145 = read paired, read reverse strand et second in pair.
# 113 = read paired, read reverse strand, mate reverse strand et first in pair
# 177 = read paired, read reverse strand, mate reverse strand et second in pair

# 113 et 177 problematique puisque aligner sur le meme brin.
# 65, 129, 81, 161, 97 et 145 insert > 3550 inclus les pairs aberentes.

quality <- subset(ecoli,ecoli$quality != '50M')
num_quality <- table(ecoli$quality)
num_quality
# 26 donc la qualite n'est pas 50M sur 8576.
# Il y a 25 non-essentiel et 1 presume essentiel.

ecoli_positif <- subset(ecoli,ecoli$length > 0 & ecoli$length < 10000)
# garde une copie par gene et enleve les longeur aberente

# calcule le chevauchement avec le gene voisin a droite
chevauchement <- c()
for (i in 1:(nrow(ecoli_positif) - 1)) {
  # Créer les séquences ens1 et ens2
  ens1 <- seq(ecoli_positif$first_pos[i], ecoli_positif$second_pos[i])
  ens2 <- seq(ecoli_positif$first_pos[i + 1], ecoli_positif$second_pos[i + 1])
  
  # Calculer le tableau de comparaison
  tab <- table(ens1 %in% ens2)
  
  if (length(tab) == 2) {
    chevauchement <- append(chevauchement,tab[2])
  }
  else {
    chevauchement <- append(chevauchement,0)
  }
}
chevauchement <- append(chevauchement,0)
ecoli_positif <- cbind(ecoli_positif, chevauchement = chevauchement)

chevauchement_gene_0 <- subset(ecoli_positif,ecoli_positif$chevauchement1 > 0)
chevauchement_gene_50 <- subset(ecoli_positif,ecoli_positif$chevauchement1 >= 50)

correspondance_l <- c()
correspondance_r <- c()
for (i in 1:nrow(ecoli_positif)) {
  print(ecoli_positif$first_pos[i])
  res <- ecoli_positif$first_pos[i] %in% ecoli_BW25113$Left.End.Position
  res2 <- ecoli_positif$second_pos[i] %in% ecoli_BW25113$Right.End.Position
  correspondance_l <- append(correspondance_l, res)
  correspondance_r <- append(correspondance_r, res2)


}
ecoli_positif <- cbind(ecoli_positif, correspondance_l = correspondance_l, correspondance_r = correspondance_r)

table(ecoli_positif$correspondance_r)
