# Fonction permettant de classer les scores par fréquence de difficultés.
modalite4 <- function(x){
  if(!is.na(x)){
    # Difficulé jamais ou très rarement répérable : {3, 4, 5}
    if(x <= 5)
      return('Très rare')
    # Difficulté peu fréquente : {6, 7, 8, 9}
    if(x >= 6 && x <= 9)
      return('Peu fréquent')
    # Difficulté fréquente : {10, 11, 12}
    if(x >= 10 && x <= 12)
      return('Fréquent')
    # Difficulté très fréquente : {13, 14, 15}
    return("Très fréquent")
  }
  return(NA)
}

fvizer <- function (sil.obj, label = FALSE, print.summary = TRUE, ...) 
{
  if (inherits(sil.obj, c("eclust", "hcut", "pam", "clara", 
                          "fanny"))) {
    df <- as.data.frame(sil.obj$silinfo$widths, stringsAsFactors = TRUE)
  }
  else if (inherits(sil.obj, "silhouette")) 
    df <- as.data.frame(sil.obj[, 1:3], stringsAsFactors = TRUE)
  else stop("Don't support an oject of class ", class(sil.obj))
  df <- df[order(df$cluster, -df$sil_width), ]
  if (!is.null(rownames(df))) 
    df$name <- factor(rownames(df), levels = rownames(df))
  else df$name <- as.factor(1:nrow(df))
  df$cluster <- as.factor(df$cluster)
  mapping <- aes_string(x = "name", y = "sil_width", color = "cluster", 
                        fill = "cluster")
  p <- ggplot(df, mapping) + geom_bar(stat = "identity") + 
    labs(y = "Silhouette width Si", x = "", title = paste0("Clusters silhouette plot ", 
                                                           "\n Average silhouette width: ",
                                                           round(mean(df$sil_width), 
                                                                 2))) + ggplot2::ylim(
                                                                   c(NA, 1)) + geom_hline(
                                                                     yintercept = mean(df$sil_width), 
                                                                     linetype = "dashed", color = "black")
  p <- ggpubr::ggpar(p, ...)
  if (!label) 
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  else if (label) 
    p <- p + theme(axis.text.x = element_text(angle = 45))
  ave <- tapply(df$sil_width, df$cluster, mean)
  n <- tapply(df$cluster, df$cluster, length)
  sil.sum <- data.frame(cluster = names(ave), size = n, ave.sil.width = round(ave, 
                                                                              2), stringsAsFactors = TRUE)
  if (print.summary) 
    print(sil.sum)
  p
}

# Transformation de la matrice noyau en matrice centrée.
kernel_centerer <- function(X){
  N_k <- nrow(X)
  one_N <- matrix(rep(1, N_k^2), nrow = N_k)/N_k
  return(X - one_N%*%X -X%*%one_N + one_N%*%X%*%one_N)
}
# Fonction permettant d'effectuer un ACP à noyau
# Noyaux : rbf, sigmoid, polynomial, cosine
k_PCA <- function(X, gamma_k=1, n_components=0, kernel_='rbf',
                  teta=0, p=2){
  if(!require(purrr))
    install.packages("purrr")
  library(purrr)
  X <- as.matrix(X)
  K <- switch (kernel_,
               # rbf kernel
               'rbf' = (- gamma_k * (X %>% stats::dist(
               ) %>% as.matrix())^2) %>% exp(),
               # sigmoid kernel
               'sigmoid' = (teta * X%*%t(X) + gamma_k) %>% tanh(),
               # polynomial kernel
               'poly' = (X%*%t(X) + gamma_k)^p,
               # cosine kernel
               'cosine' = (X%*%t(X)) %>% cov()
  ) %>% kernel_centerer() %>% eigen(symmetric = TRUE)
  # récupération des indices des valeurs propres non nulles
  non_zeros_ind <- which(K$values != 0)
  # récupération des vecteurs propres dont les valeurs propre
  # sont non nulles
  vects <- K$vectors[, non_zeros_ind]
  n_ <- ncol(X)
  # Labelisation des composantes principales
  colnames(vects) <- paste('CP', 1:nrow(X), sep = '')
  if(n_components == 0)
    return(vects[, 1:n_])
  return(vects[, 1:min(n_, n_components)])
}

# *****************Fin de la définition des fonction **********


# Chargement des données
donnees <- read.csv("~/Desktop/M2-SAAD/stage_proj/data.csv", sep = ",",
                    na.strings = c(NA, ""), row.names = 1)
# Nommons les variables
colnames(donnees) <- c("NbreEfts", "GenreEft", "AgeEft",              
                       "RangEft", "AllaitSeinMat", "AllaitLaitMat",       
                       "AllaitLaitMatHypoall", "AllaitLaitEpai", "AllaitLaitCrois",    
                       "LaitVache", "DureeAllaitSeinMat", "EftAccueilCol",       
                       "NbreRepasCol", "ManquePlaisir", "Neophobie",           
                       "Selectivite", "ManqueAppetit", "Preference",          
                       "Explication", "Coercision", "Contingence",         
                       "BiaisPercepParent", "AnxieteParent", "NaisTermeGros",       
                       "NaisCmbienSem", "EftHospitalPlus1Nuit", "IMCEstime",           
                       "TailleAge", "PoidsAge", "PoidsTaille",         
                       "EftAll", "NbreDmgeDep8Nais", "VieCouple",           
                       "AgeMere", "AgePere", "NivoEtuMere",     
                       "MereEtu", "MereSalariee", "CSPMere", "NivoEtuPere",     
                       "PereEtu", "PereSalarie", "CSPPere", "RevenuAnFoyer")
# attach(donnees)
# detach(donnees)
# Label encoding
# Conversion du genre des enfant
# 0 pour masculin et 1 féminin
donnees$RangEft <- factor(donnees$RangEft)

# Traitement des données manquantes
library(missMDA)
res.X_imputes <- imputeFAMD(donnees, ncp = 3)
X_impute <- res.X_imputes$completeObs
attach(X_impute)

# Extraction des données numériques des enfants: colonnes
donnees_enfants <- X_impute[, c(1, 3, 11, 13:23, 27:30, 34, 35)]
# donnees_enfants <- donnees_enfants[-c(91, 263, 337:339), ]
# donnees_enfants <- donnees_enfants[, -c(9:14)]
# Number of observation
n <- nrow(donnees_enfants)
acp_2_comp <- k_PCA(donnees_enfants, gamma_k = 1/n,
                    n_components = 2, kernel_ = 'rbf')
# Nombre d'observations
n <- nrow(acp_2_comp)
# How many centers would you take to our kmeans clustering?
# Total sum of square wihtin each cluster.
# Calcul du nombre idéal pour effectuer le K-means

# Calcul de
totwithinssList <- NULL
for (i in 1:10) {
  centers = acp_2_comp[sample(1:n, i, replace = FALSE), ]
  if(i==1)
    centers = data.frame(centers[1], centers[2])
  km <- kmeans(acp_2_comp, centers = centers, algorithm = 'Hartigan-Wong',
               iter.max = 300)
  totwithinssList[i] <- km$tot.withinss
}

# Graphique de la méthode de elbow.
plot(totwithinssList, type = 'o', col = '#626567', lty = 1,
     ylab = "Distorsions", xlab = "Number of clusters (k)", lwd=2,
     main = "")
points(totwithinssList, pch=16)
#
abline(v=4, lty=3, col='#626567', lwd=3)
grid()
legend('topright', col = c(1, '#626567'), 
       pch = c(19, NA), lty = c(NA, 3),
       legend = c('Distorsions', 'Optimal k'), bg = '#f4f6f6',
       box.lwd = NA, cex = 1, lwd = 3)

# Kmeans avec 04 clusters
library(clustertend)
# library(ggplot2)
library(factoextra)
# Calcul de la métrique de Hopkins
# hopkins(donnees_enfants, n=nrow(donnees_enfants)-1)

# Retrait des individus 91, 263, 337:339
# acp_2_comp <- acp_2_comp[-c(91, 263, 337:339), ]

# Vecteur des couleurs
# colors_ = c("#149414", '#955628', "#6C0277", "#007FFF")
colors_ <- c('#e5e7e9', '#bdc3c7', '#909497', '#626567')
# Vecteur des puces
puces = 15:18
# Nombre de clusters
cluster_number <- 4
km_4_clusters <- eclust(acp_2_comp, "kmeans", k = cluster_number, 
                        nstart = 25, graph = F)
# Plot de silhouette des Clusters
fvizer(km_4_clusters, ggtheme=theme_classic(),
                palette=colors_[1:cluster_number]) 

# Affichage des clusters
plot(acp_2_comp[km_4_clusters$cluster==1, ],
     pch = puces[1], col = colors_[1], lwd = 2, 
     ylab = "CP2", xlab = "CP1",
     xlim = c(min(c(acp_2_comp[, 1], km_4_clusters$centers[, 1]))-.01, 
              max(c(acp_2_comp[, 1], km_4_clusters$centers[, 1]))+.02),
     ylim = c(min(c(acp_2_comp[, 2], km_4_clusters$centers[, 2]))-.01, 
              max(c(acp_2_comp[, 2], km_4_clusters$centers[, 2]))+.01)
     )
nameLeng <- c("Cluster 1", rep("0", cluster_number-1))
for (c in 2:cluster_number) {
  points(acp_2_comp[km_4_clusters$cluster==c, ],
         pch = puces[c], col = colors_[c], lwd=2)
  nameLeng[c] <- paste("Cluster", c, sep = " ")
}

# # Affichage des individus à retirer
# points(acp_2_comp[c(91, 263, 337:339), ], pch = 19, col = 1,
#        lwd = 4)
# text(acp_2_comp[c(91, 263, 337:339), ] - c(rep(0.006, 3), -0.006, 0.006),
#      labels = factor(c(91, 263, 337:339)),
#      col = 1)

# Les centres des Clusters
points(km_4_clusters$centers, pch = 9, lwd=3, col = "#1b2631")
# color="Gris perle" grid
grid(col = '#CECECE')
# La légende du nuage de points
legend('topright', col = c( '#bdc3c7', '#909497', '#626567', '#e5e7e9', "#1b2631"), 
       pch = c(16, 17, 18, 15, 9),
       legend = c(nameLeng[1:cluster_number], "Centers"), bg = '#f2f3f4',
       box.lwd = NA, cex = .8)

# Laison des individus de chaque cluster à leur centre
for (c_ in 1:cluster_number) {
  segments(acp_2_comp[km_4_clusters$cluster==c_, 1],
           acp_2_comp[km_4_clusters$cluster==c_, 2],
           km_4_clusters$centers[c_, 1],
           km_4_clusters$centers[c_, 2],
           col = colors_[c_], lwd=1, lty = 1)
}

abline(h=0, lty=2, lwd=2); abline(v=0, lty=2, lwd=2)
# Fin plotting

# ***************************************************
# Détermination des caractères les plus importants **
# et des parangons de chaque cluster                *
# ***************************************************
# Centre d'inertie du nuage de points

parangX_impo <- function(X, clusterin){
  X_bar_j <- apply(X = X, MARGIN = 2, FUN = mean)
  s_j_2 <- apply(X = X, MARGIN = 2, FUN = var)
  n_d <- nrow(X)
  pvalue_jg <- NULL
  parangons <- NULL
  for(g in 1:clusterin$nbclust){
    cl <- paste0('clus', g)
    # Les individus du cluster g
    cls <- X[clusterin$cluster==g, ]
    # Taille du cluste g.
    n_g <- nrow(cls)
    # Le centre d'inertie du cluster g.
    X_bar_jg <- apply(X = cls, MARGIN = 2, FUN = mean)
    # Détermination du parangon du cluster g.
    # index du parangons dans les données.
    pg <- (rbind(X_bar_jg, cls) %>% dist() %>% as.matrix(
    ))[1, -1] %>% which.min() %>%names()
    parangons[[paste0('prgon', g, '_is_', pg)]] <- X[pg, ]
    z_j_g <- (X_bar_jg - X_bar_j)/(s_j_2*(n_d - n_g)/(n_g*(n_d - 1)))^0.5
    # P-valeur donnant l'importance de la variable j dans le cluster g
    pvalue_jg[[cl]] <- sort(2*(1 - pnorm(abs(z_j_g))))
  }
  return(list('parangons'=parangons, 'p_value'=pvalue_jg))
}

pX_imp <- parangX_impo(donnees_enfants, km_4_clusters)
pX_imp$parangons
pX_imp$p_value

# # ****************** Min, Max, Mean et SD de chaque cluster********
# X_clus <- NULL
# for(clus in 1:4){
#   X_pas <- NULL
#   for(fn in c('min', 'max', 'mean', 'sd')){
#     X_pas[[fn]] <- round(apply(donnees_enfants[km_4_clusters$cluster==clus, ],
#                          MARGIN = 2, FUN = fn), 2)
#   }
#   X_clus[[paste0('cluster', clus)]] <- as.data.frame(X_pas)
# }
# X_clus
# 
# # Correlation between numeric variables
# cor_num_data <- cor(donnees_enfants)
# cor_vd <- cor(donnees_enfants[, 5:8], method = 'spearman')
# library(ggplot2)
# library('corrplot')
# 
# corrplot(cor_num_data, type = 'upper', method = 'square')
# corrplot(cor_vd, type = 'upper', method = 'number')
# 
# # ***************************************************
# # **************** LES ARBRES DE DÉCISION ***********
# # ***************************************************
# library(rpart)
# library(rpart.plot)
# 
# AppetitMod <- NULL
# PlaisirMod <- NULL
# NeophobieMod <- NULL
# SelectiviteMod <- NULL
# PreferenceMod <- NULL
# ExplicationMod <- NULL
# CoercisionMod <- NULL
# ContingenceMod <- NULL
# 
# X_impute1 <- X_impute[-c(91, 263, 337:339), ]
# # X_impute1<- X_impute[-sortants, ]
# n <- nrow(X_impute1)
# 
# for (i in 1:n){
#   AppetitMod[i] <- modalite4(X_impute1$ManqueAppetit[i])
#   PlaisirMod[i] <- modalite4(X_impute1$ManquePlaisir[i])
#   NeophobieMod[i] <- modalite4(X_impute1$Neophobie[i])
#   SelectiviteMod[i] <- modalite4(X_impute1$Selectivite[i])
#   PreferenceMod[i] <- modalite4(X_impute1$Preference[i])
#   ExplicationMod[i] <- modalite4(X_impute1$Explication[i])
#   CoercisionMod[i] <- modalite4(X_impute1$Coercision[i])
#   ContingenceMod[i] <- modalite4(X_impute1$Contingence[i])
# }
# 
# vdData <- X_impute1[, 14:17]
# d1 <- cbind(vdData[, 1:3], AppetitMod)
# d2 <- cbind(vdData[, 2:4], PlaisirMod)
# d3 <- cbind(vdData[, c(1,3,4)], NeophobieMod)
# d4 <- cbind(vdData[, c(1,2,4)], SelectiviteMod)
# 
# # ***************************
# #** Variable ManqueAppetit **
# # ***************************
# arbre1 <- rpart(AppetitMod ~ ., data = d1)
# rpart.plot(arbre1, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
#            # Fréquent, Peu fréquent, Très fréquent, Très rare.
#            # '#A4A4A4','#D8D8D8', '#6E6E6E', '#F2F2F2' Gris
#            box.palette = list('#FEA347','#FEE347', '#FF4901', '#1FA055'),
#            shadow.col = 'grey')
# rpart.rules(arbre1)
# # COMMENTAIRE : Manque appetit avec les 4 VD
# # Un enfant a très fréquemment un manque d'appetit lorsqu'il a souvent un manque
# # de plaisir avec un score supérieur ou égal à 11. Et ils font 7% de l'échantillon.
# # Il peut être considéré avoir un manque d'appétit fréquent, lorsqu'il a un manque
# # de plaisir peu fréquent ou fréquent (8<= ManquePlaisir <=11) ou il est très
# # sélectif (Selectivite >=13) et a rarement ou peu fréquemment un manque de plaisir
# # (4<= ManquePlaisir <= 8). Et parmi les enfants étudiés, 22% ont un manque d'appétit.
# # Nous pouvons déduire aussi qu'à partir d'un score de manque de plaisir supérieur
# # ou égal à 8, qu'un enfant a un manque d'appétit fréquent.
# # Une forte partie de l'échantillon (64%) a peu fréquemment un manque d'appétit.
# # Ce manque d'appétit peu fréquent peut se décider chez un enfant, lorsqu'il ne
# # manque pas très souvent de plaisir (ManquePlaisir < 4) et est au moins peu
# # sélectif (Selectivite >=8) ou bien il a un manque de plaisir peu fréquent ou très
# # rare (4<=ManquePlaisir<8) et n'est pas très sélectif (Selectivite <13).
# # Enfin un enfant a très fréquemment d'appétit s'il a un manque de plaisir très
# # rare (ManquePlaisir < 4) et est au plus peu sélectif (Selectivite < 8). Et ils sont
# # peu nombreux (6%) dans la population étudiée.
# 
# # Taux d'erreur de l'arbre de classification de la variable ManqueAppetit
# mean(residuals(arbre1))
# 
# plotcp(arbre1)
# arbre11 <- prune(arbre1, cp = 0.028)
# rpart.plot(arbre11, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
#            # Fréquent, Peu fréquent, Très fréquent, Très rare.
#            # box.palette = list('#A4A4A4','#D8D8D8', '#6E6E6E', '#F2F2F2'),
#            shadow.col = 'grey')
# rpart.rules(arbre11)
# mean(residuals(arbre11))
# 
# # ***************************
# #** Variable ManquePlaisir **
# # ***************************
# arbre2 <- rpart(PlaisirMod ~ ., data = d2)
# rpart.plot(arbre2, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
#            # Fréquent, Peu fréquent, Très fréquent, Très rare.
#            # '#A4A4A4','#D8D8D8', '#6E6E6E', '#F2F2F2'
#            box.palette = list('#FEA347','#FEE347', #'#FF4901',
#                               '#1FA055'),
#            shadow.col = 'grey')
# rpart.rules(arbre2)
# # Un manque de plaisir fréquent est identifé chez un enfant lorsqu'il a un manque
# # d'appétit fréquent ou très fréquent (ManqueAppetit >=11) et est souvent ou
# # très souvent sélectif (9<= Selectivite).
# 
# # La majorité de (60%) des enfants ont peu fréquemment un manque de plaisir
# # lorsqu'ils ont un manque d'appétit très fréquent (ManqueAppetit >=11) ou ils n'ont
# # pas souvent ou peu souvent d'appétit et ne sont pas très frquemment
# # néophobes.
# 
# # On remarque 33% des individus prennent très frequemment du plaisir à manger
# # lorsqu'ils ont souvent un appétit (ManqueAppetit < 7) ils représentent
# #  25% des 33% ou ils ont un peu ou fréquemment d'appétit (7< ManqueAppetit <11)
# # et sont très rarement néophobes soit ils sont souvent néophobes (Neophobie >=11).
# 
# # Taux d'erreur de l'arbre de classification de la variable ManqueAppetit
# mean(residuals(arbre2))
# 
# plotcp(arbre2)
# arbre22 <- prune(arbre2, cp = 0.026)
# rpart.plot(arbre22, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
#            # Fréquent, Peu fréquent, Très fréquent, Très rare.
#            # box.palette = list('#A4A4A4','#D8D8D8', '#6E6E6E', '#F2F2F2'),
#            shadow.col = 'grey')
# rpart.rules(arbre22)
# mean(residuals(arbre22))
# 
# # ***************************
# #** Variable Neophobie **
# # ***************************
# arbre3 <- rpart(NeophobieMod ~ ., data = d3)
# rpart.plot(arbre3, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
#            # Fréquent, Peu fréquent, Très fréquent, Très rare.
#            box.palette = list('#FEA347','#FEE347', '#FF4901', '#1FA055'),
#            shadow.col = 'grey')
# rpart.rules(arbre3)
# # Les enfants sont très fréquemment néophobes lorsqu'ils sont très souvent sélectifs
# # (Selectivite >=13) et ils constituent 8% de l'échantillon. A l'inverse, les enfants
# # sont très rarement néophobes s'ils sont très rarement sélectifs (Selectivite < 6).
# # Ils font 17% de la population étudiée.
# # Certains d'entre eux (20%), fréquemment néophobes lorsqu'ils sont aussi
# # fréquemment sélectifs. Et le reste de ces enfants (55%) qui est la grande
# # partie de l'échantillon sont peu fréquemment nophobes lorsqu'ils sont
# # peu fréquemment sélectifs uniquement ou ils manquent très rarement ou
# #  peu fréquemment de plaisir à manger (ManquePlaisir < 10) et sont peu
# # souvent sélectifs (Sélectivite=9).
# 
# # Taux d'erreur de l'arbre de classification de la variable ManqueAppetit
# mean(residuals(arbre3))
# 
# plotcp(arbre3)
# arbre33 <- prune(arbre3, cp = 0.044)
# rpart.plot(arbre33, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
#            # Fréquent, Peu fréquent, Très fréquent, Très rare.
#            # box.palette = list('#A4A4A4','#D8D8D8', '#6E6E6E', '#F2F2F2'),
#            shadow.col = 'grey')
# rpart.rules(arbre33)
# mean(residuals(arbre33))
# 
# # ***************************
# #** Variable Selectivite **
# # ***************************
# arbre4 <- rpart(SelectiviteMod ~ ., data = d4)
# rpart.plot(arbre4, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
#            # Fréquent, Peu fréquent, Très fréquent, Très rare.
#            box.palette = list('#FEA347','#FEE347', '#FF4901', '#1FA055'),
#            shadow.col = 'grey')
# rpart.rules(arbre4)
# # Taux d'erreur de l'arbre de classification de la variable ManqueAppetit
# mean(residuals(arbre4))
# 
# plotcp(arbre4)
# arbre44 <- prune(arbre4, cp = 0.033)
# rpart.plot(arbre44, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
#            # Fréquent, Peu fréquent, Très fréquent, Très rare.
#            # box.palette = list('#A4A4A4','#D8D8D8', '#6E6E6E', '#F2F2F2'),
#            shadow.col = 'grey')
# rpart.rules(arbre44)
# mean(residuals(arbre44))
# 
# # ****************************************************
# # ****************************************************
# #  Les arbres de classification sans les stratégies***
# #  parentales                                      ***
# # ****************************************************
# # ****************************************************
# # c(2:4, 11, 13, 24, 26, 30:36, -33, 38, 40, 42, 44)
# # -c(14:23, 39, 43, 27:30)
# vdDataSansStratP <- X_impute1[, c(2:4, 11, 13, 24, 26, 30:32,
#                                     34:36, 38, 40, 42, 44)]
# # 14:ManquePlaisir 15:Neophobie 16:Selectivite 17:ManqueAppetit
# dd1 <- cbind(vdDataSansStratP, AppetitMod)
# dd2 <- cbind(vdDataSansStratP, PlaisirMod)
# dd3 <- cbind(vdDataSansStratP, NeophobieMod)
# dd4 <- cbind(vdDataSansStratP, SelectiviteMod)
# 
# # ***************************
# #** Variable ManqueAppetit **
# # ***************************
# arbre1 <- rpart(AppetitMod ~ ., data = dd1, maxdepth=5)
# rpart.plot(arbre1, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
#            # Fréquent, Peu fréquent, Très fréquent, Très rare.
#            #box.palette = list('#A4A4A4','#D8D8D8', '#6E6E6E', '#F2F2F2'),
#            shadow.col = 'grey')
# rpart.rules(arbre1)
# mean(residuals(arbre1))
# 
# # plotcp(arbre1)
# # arbre11 <- prune(arbre1, cp = 0.027)
# # rpart.plot(arbre11, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
# #            # Fréquent, Peu fréquent, Très fréquent, Très rare.
# #            # box.palette = list('#A4A4A4','#D8D8D8', '#6E6E6E', '#F2F2F2'),
# #            shadow.col = 'grey')
# # rpart.rules(arbre11)
# # mean(residuals(arbre11))
# 
# # ***************************
# #** Variable ManquePlaisir **
# # ***************************
# arbre2 <- rpart(PlaisirMod ~ ., data = dd2, maxdepth=3)
# # rpart.plot(arbre2, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
# #            # Fréquent, Peu fréquent, Très fréquent, Très rare.
# #            #box.palette = list('#A4A4A4','#D8D8D8', '#6E6E6E', '#F2F2F2'),
# #            shadow.col = 'grey')
# rpart.rules(arbre2)
# # Taux d'erreur de l'arbre de classification de la variable ManqueAppetit
# mean(residuals(arbre2))
# 
# # plotcp(arbre2)
# # arbre22 <- prune(arbre2, cp = 0.023)
# # rpart.plot(arbre22, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
# #            # Fréquent, Peu fréquent, Très fréquent, Très rare.
# #            # box.palette = list('#A4A4A4','#D8D8D8', '#6E6E6E', '#F2F2F2'),
# #            shadow.col = 'grey')
# # rpart.rules(arbre22)
# # mean(residuals(arbre22))
# 
# # ***************************
# #** Variable Neophobie **
# # ***************************
# arbre3 <- rpart(NeophobieMod ~ ., data = dd3, maxdepth=3)
# # rpart.plot(arbre3, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
# #            # Fréquent, Peu fréquent, Très fréquent, Très rare.
# #            #box.palette = list('#A4A4A4','#D8D8D8', '#6E6E6E', '#F2F2F2'),
# #            shadow.col = 'grey', cex = .3)
# rpart.rules(arbre3)
# # Taux d'erreur de l'arbre de classification de la variable ManqueAppetit
# mean(residuals(arbre3))
# 
# # plotcp(arbre3)
# # arbre33 <- prune(arbre3, cp = 0.066)
# # rpart.plot(arbre33, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
# #            # Fréquent, Peu fréquent, Très fréquent, Très rare.
# #            # box.palette = list('#A4A4A4','#D8D8D8', '#6E6E6E', '#F2F2F2'),
# #            shadow.col = 'grey')
# # rpart.rules(arbre33)
# # mean(residuals(arbre33))
# 
# # ***************************
# #** Variable Selectivite **
# # ***************************
# arbre4 <- rpart(SelectiviteMod ~ ., data = dd4, maxdepth=4)
# # rpart.plot(arbre4, legend.x = NA, extra = 100, type = 1, branch.lty = 3, nn = TRUE,
# #            # Fréquent, Peu fréquent, Très fréquent, Très rare.
# #            #box.palette = list('#A4A4A4','#D8D8D8', '#6E6E6E', '#F2F2F2'),
# #            shadow.col = 'grey')
# rpart.rules(arbre4)
# # Taux d'erreur de l'arbre de classification de la variable ManqueAppetit
# mean(residuals(arbre4))
# 
# # plotcp(arbre4)
# # arbre44 <- prune(arbre4, cp = 0.036)
# # rpart.plot(arbre44, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
# #            # Fréquent, Peu fréquent, Très fréquent, Très rare.
# #            # box.palette = list('#A4A4A4','#D8D8D8', '#6E6E6E', '#F2F2F2'),
# #            shadow.col = 'grey')
# # rpart.rules(arbre44)
# # mean(residuals(arbre44))
# 
# # ********************************************************************
# # ********************************************************************
# # RELATIONS ENTRE DIFFICULTES ALIMENTAIRES ET STRATEGIES PARENTALES **
# # ********************************************************************
# # ********************************************************************
# vpData <- X_impute[-c(91, 263, 337:339), c(18:21)]
# vpData <- round(vpData)
# ddd1 <- cbind(vpData, AppetitMod)
# ddd2 <- cbind(vpData, PlaisirMod)
# ddd3 <- cbind(vpData, NeophobieMod)
# ddd4 <- cbind(vpData, SelectiviteMod)
# 
# # ***************************
# #** Variable ManqueAppetit **
# # ***************************
# arbre1 <- rpart(AppetitMod ~ ., data = ddd1)
# rpart.plot(arbre1, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
#            # Fréquent, Peu fréquent, Très fréquent, Très rare.
#            box.palette = list('#FEA347','#FEE347', '#FF4901', '#1FA055'),
#            shadow.col = 'grey')
# rpart.rules(arbre1)
# mean(residuals(arbre1))
# 
# # plotcp(arbre1)
# # arbre11 <- prune(arbre1, cp = 0.027)
# # rpart.plot(arbre11, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
# #            # Fréquent, Peu fréquent, Très fréquent, Très rare.
# #            # box.palette = list('#A4A4A4','#D8D8D8', '#6E6E6E', '#F2F2F2'),
# #            shadow.col = 'grey')
# # rpart.rules(arbre11)
# # mean(residuals(arbre11))
# 
# # ***************************
# #** Variable ManquePlaisir **
# # ***************************
# arbre2 <- rpart(PlaisirMod ~ ., data = ddd2)
# rpart.plot(arbre2, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
#            # Fréquent, Peu fréquent, Très fréquent, Très rare.
#            box.palette = list(#'#FEA347',
#                               '#FEE347', #'#FF4901',
#                               '#1FA055'),
#            shadow.col = 'grey')
# rpart.rules(arbre2)
# # Taux d'erreur de l'arbre de classification de la variable ManqueAppetit
# mean(residuals(arbre2))
# 
# # plotcp(arbre2)
# # arbre22 <- prune(arbre2, cp = 0.023)
# # rpart.plot(arbre22, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
# #            # Fréquent, Peu fréquent, Très fréquent, Très rare.
# #            # box.palette = list('#A4A4A4','#D8D8D8', '#6E6E6E', '#F2F2F2'),
# #            shadow.col = 'grey')
# # rpart.rules(arbre22)
# # mean(residuals(arbre22))
# 
# # ***************************
# #** Variable Neophobie **
# # ***************************
# arbre3 <- rpart(NeophobieMod ~ ., data = ddd3)
# rpart.plot(arbre3, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
#            # Fréquent, Peu fréquent, Très fréquent, Très rare.
#            box.palette = list('#FEA347','#FEE347', '#FF4901', '#1FA055'),
#            shadow.col = 'grey')
# rpart.rules(arbre3)
# # Taux d'erreur de l'arbre de classification de la variable ManqueAppetit
# mean(residuals(arbre3))
# 
# # plotcp(arbre3)
# # arbre33 <- prune(arbre3, cp = 0.066)
# # rpart.plot(arbre33, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
# #            # Fréquent, Peu fréquent, Très fréquent, Très rare.
# #            # box.palette = list('#A4A4A4','#D8D8D8', '#6E6E6E', '#F2F2F2'),
# #            shadow.col = 'grey')
# # rpart.rules(arbre33)
# # mean(residuals(arbre33))
# 
# # ***************************
# #** Variable Selectivite **
# # ***************************
# arbre4 <- rpart(SelectiviteMod ~ ., data = ddd4)
# rpart.plot(arbre4, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
#            # Fréquent, Peu fréquent, Très fréquent, Très rare.
#            box.palette = list('#FEA347','#FEE347'#, '#FF4901', '#1FA055'
#                               ),
#            shadow.col = 'grey', cex = 0.5)
# rpart.rules(arbre4)
# # Taux d'erreur de l'arbre de classification de la variable ManqueAppetit
# mean(residuals(arbre4))
# 
# # plotcp(arbre4)
# # arbre44 <- prune(arbre4, cp = 0.036)
# # rpart.plot(arbre44, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
# #            # Fréquent, Peu fréquent, Très fréquent, Très rare.
# #            # box.palette = list('#A4A4A4','#D8D8D8', '#6E6E6E', '#F2F2F2'),
# #            shadow.col = 'grey')
# # rpart.rules(arbre44)
# # mean(residuals(arbre44))
# 
# # RELATIONS ENTRE STRATEGIES PARENTALES et DIFFICULTES PARENTALES
# dddd1 <- cbind(PreferenceMod, vdData)
# dddd2 <- cbind(ExplicationMod, vdData)
# dddd3 <- cbind(CoercisionMod, vdData)
# dddd4 <- cbind(ContingenceMod, vdData)
# 
# # ***************************
# #** Variable Preference **
# # ***************************
# arbre1 <- rpart(PreferenceMod ~ ., data = dddd1)
# rpart.plot(arbre1, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
#            # Fréquent, Peu fréquent, Très fréquent, Très rare.
#            box.palette = list('#FEA347','#FEE347'#, '#FF4901', '#1FA055'
#                               ),
#            shadow.col = 'grey')
# rpart.rules(arbre1)
# mean(residuals(arbre1))
# 
# # plotcp(arbre1)
# # arbre11 <- prune(arbre1, cp = 0.027)
# # rpart.plot(arbre11, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
# #            # Fréquent, Peu fréquent, Très fréquent, Très rare.
# #            # box.palette = list('#A4A4A4','#D8D8D8', '#6E6E6E', '#F2F2F2'),
# #            shadow.col = 'grey')
# # rpart.rules(arbre11)
# # mean(residuals(arbre11))
# 
# # ***************************
# #** Variable Explication **
# # ***************************
# arbre2 <- rpart(ExplicationMod ~ ., data = dddd2)
# rpart.plot(arbre2, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
#            # Fréquent, Peu fréquent, Très fréquent, Très rare.
#            box.palette = list('#FEA347','#FEE347', '#FF4901'#, '#1FA055'
#            ),
#            shadow.col = 'grey')
# rpart.rules(arbre2)
# # Taux d'erreur de l'arbre de classification de la variable ManqueAppetit
# mean(residuals(arbre2))
# 
# # plotcp(arbre2)
# # arbre22 <- prune(arbre2, cp = 0.023)
# # rpart.plot(arbre22, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
# #            # Fréquent, Peu fréquent, Très fréquent, Très rare.
# #            # box.palette = list('#A4A4A4','#D8D8D8', '#6E6E6E', '#F2F2F2'),
# #            shadow.col = 'grey')
# # rpart.rules(arbre22)
# # mean(residuals(arbre22))
# 
# # ***************************
# #** Variable Coercision **
# # ***************************
# arbre3 <- rpart(CoercisionMod ~ ., data = dddd3)
# rpart.plot(arbre3, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
#            # Fréquent, Peu fréquent, Très fréquent, Très rare.
#            box.palette = list('#FEA347','#FEE347', #'#FF4901',
#                               '#1FA055'),
#            shadow.col = 'grey')
# rpart.rules(arbre3)
# # Taux d'erreur de l'arbre de classification de la variable ManqueAppetit
# mean(residuals(arbre3))
# 
# # plotcp(arbre3)
# # arbre33 <- prune(arbre3, cp = 0.066)
# # rpart.plot(arbre33, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
# #            # Fréquent, Peu fréquent, Très fréquent, Très rare.
# #            # box.palette = list('#A4A4A4','#D8D8D8', '#6E6E6E', '#F2F2F2'),
# #            shadow.col = 'grey')
# # rpart.rules(arbre33)
# # mean(residuals(arbre33))
# 
# # ***************************
# #** Variable Contingence **
# # ***************************
# arbre4 <- rpart(ContingenceMod ~ ., data = dddd4)
# rpart.plot(arbre4, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
#            # Fréquent, Peu fréquent, Très fréquent, Très rare.
#            box.palette = list(#'#FEA347',
#                               '#FEE347', #'#FF4901', 
#                               '#1FA055'),
#            shadow.col = 'grey')
# rpart.rules(arbre4)
# # Taux d'erreur de l'arbre de classification de la variable ManqueAppetit
# mean(residuals(arbre4))
# 
# # plotcp(arbre4)
# # arbre44 <- prune(arbre4, cp = 0.036)
# # rpart.plot(arbre44, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
# #            # Fréquent, Peu fréquent, Très fréquent, Très rare.
# #            # box.palette = list('#A4A4A4','#D8D8D8', '#6E6E6E', '#F2F2F2'),
# #            shadow.col = 'grey')
# # rpart.rules(arbre44)
# # mean(residuals(arbre44))
# 
# # *******************************
# # ENTRE LES STRATEGIES PARENTALES
# # *******************************
# ddddd1 <- cbind(vpData[, -1], PreferenceMod)
# ddddd2 <- cbind(vpData[, -2], ExplicationMod)
# ddddd3 <- cbind(vpData[, -3], CoercisionMod)
# ddddd4 <- cbind(vpData[, -4], ContingenceMod)
# 
# # ***************************
# #** Variable Preference **
# # ***************************
# arbre1 <- rpart(PreferenceMod ~ ., data = ddddd1)
# rpart.plot(arbre1, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
#            # Fréquent, Peu fréquent, Très fréquent, Très rare.
#            box.palette = list('#FEA347','#FEE347'#, '#FF4901', '#1FA055'
#            ),
#            shadow.col = 'grey')
# rpart.rules(arbre1)
# mean(residuals(arbre1))
# 
# # plotcp(arbre1)
# # arbre11 <- prune(arbre1, cp = 0.027)
# # rpart.plot(arbre11, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
# #            # Fréquent, Peu fréquent, Très fréquent, Très rare.
# #            # box.palette = list('#A4A4A4','#D8D8D8', '#6E6E6E', '#F2F2F2'),
# #            shadow.col = 'grey')
# # rpart.rules(arbre11)
# # mean(residuals(arbre11))
# 
# # ***************************
# #** Variable Explication **
# # ***************************
# arbre2 <- rpart(ExplicationMod ~ ., data = ddddd2)
# rpart.plot(arbre2, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
#            # Fréquent, Peu fréquent, Très fréquent, Très rare.
#            box.palette = list('#FEA347','#FEE347', '#FF4901'#, '#1FA055'
#                               ),
#            shadow.col = 'grey')
# rpart.rules(arbre2)
# # Taux d'erreur de l'arbre de classification de la variable ManqueAppetit
# mean(residuals(arbre2))
# 
# # plotcp(arbre2)
# # arbre22 <- prune(arbre2, cp = 0.023)
# # rpart.plot(arbre22, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
# #            # Fréquent, Peu fréquent, Très fréquent, Très rare.
# #            # box.palette = list('#A4A4A4','#D8D8D8', '#6E6E6E', '#F2F2F2'),
# #            shadow.col = 'grey')
# # rpart.rules(arbre22)
# # mean(residuals(arbre22))
# 
# # ***************************
# #** Variable Coercision **
# # ***************************
# arbre3 <- rpart(CoercisionMod ~ ., data = ddddd3)
# rpart.plot(arbre3, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
#            # Fréquent, Peu fréquent, Très fréquent, Très rare.
#            box.palette = list('#FEA347','#FEE347', #'#FF4901',
#                               '#1FA055'),
#            shadow.col = 'grey')
# rpart.rules(arbre3)
# # Taux d'erreur de l'arbre de classification de la variable ManqueAppetit
# mean(residuals(arbre3))
# 
# # plotcp(arbre3)
# # arbre33 <- prune(arbre3, cp = 0.066)
# # rpart.plot(arbre33, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
# #            # Fréquent, Peu fréquent, Très fréquent, Très rare.
# #            # box.palette = list('#A4A4A4','#D8D8D8', '#6E6E6E', '#F2F2F2'),
# #            shadow.col = 'grey')
# # rpart.rules(arbre33)
# # mean(residuals(arbre33))
# 
# # ***************************
# #** Variable Contingence **
# # ***************************
# arbre4 <- rpart(ContingenceMod ~ ., data = ddddd4)
# rpart.plot(arbre4, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
#            # Fréquent, Peu fréquent, Très fréquent, Très rare.
#            box.palette = list(#'#FEA347',
#              '#FEE347', #'#FF4901', 
#              '#1FA055'),
#            shadow.col = 'grey')
# rpart.rules(arbre4)
# # Taux d'erreur de l'arbre de classification de la variable ManqueAppetit
# mean(residuals(arbre4))
# 
# # plotcp(arbre4)
# # arbre44 <- prune(arbre4, cp = 0.036)
# # rpart.plot(arbre44, legend.x = NA, extra = 100, type = 4, branch.lty = 3, nn = TRUE,
# #            # Fréquent, Peu fréquent, Très fréquent, Très rare.
# #            # box.palette = list('#A4A4A4','#D8D8D8', '#6E6E6E', '#F2F2F2'),
# #            shadow.col = 'grey')
# # rpart.rules(arbre44)
# # mean(residuals(arbre44))
# 
# # *******************************************************
# # ****************** REGRESSION MULTINOMIALE ************
# # *******************************************************
# X_impute1 <- X_impute[-c(91, 263, 337:339), ]
# X_impute1 <- subset(X_impute1, select = -c(CSPMere, CSPPere))
# 
# # Régression pour ManqueAppetit
# rDataManqueApp <- cbind(AppetitMod, X_impute1)
# rDataManqueApp <- subset(rDataManqueApp, select = -c(ManqueAppetit))
# 
# # mise en oeuvre de la regression
# library(MASS)
# library(car)
# rgAppetit <- polr(AppetitMod~., data = rDataManqueApp, method = 'logistic')
# summary(rgAppetit)
# 
# 
# # Sélection de variables
# library(stats)
# rgAppetit_sel <- step(rgAppetit, k = 2, direction = 'both')
# # Comparaison de modèle avant et après sélection de variables
# anova(rgAppetit, rgAppetit_sel)
# # A la suite de la comparaison de modèle avec la fonction anova,
# # nous avons une p-valeur = 0.1232089 ce qui indique qu'il n'y a
# # pas de différence entre les deux modèles.
# summary(rgAppetit_sel)
# zAppetit_sel = summary(rgAppetit_sel)$coeff / summary(rgAppetit_sel)$standard.errors
# pval_Appetit_sel = 2 * (1 - pnorm(abs(zAppetit_sel), 0, 1))
# pval_Appetit_sel
# prAppet <- predict(rgAppetit_sel)
# mcAppet <- table(AppetitMod, prAppet)
# tauxAppet <- 1 - sum(diag(mcAppet))/sum(mcAppet)
# tauxAppet
# 
# # **************************************************
# # Régression pour ManquePlaisir
# rDataManquePlaisir <- cbind(PlaisirMod, X_impute1)
# rDataManquePlaisir <- subset(rDataManquePlaisir, select = -c(ManquePlaisir))
# 
# # mise en oeuvre de la regression
# library(nnet)
# rgPlaisir <- multinom(PlaisirMod~., data = rDataManquePlaisir)
# summary(rgPlaisir)
# zPlaisir = summary(rgPlaisir)$coeff / summary(rgPlaisir)$standard.errors
# pval_Plaisir = 2 * (1 - pnorm(abs(zPlaisir), 0, 1))
# pval_Plaisir
# 
# # Sélection de variables
# library(stats)
# rgPlaisir_sel <- step(rgPlaisir, k = 2, direction = 'both')
# # Comparaison de modèle avant et après sélection de variables
# anova(rgPlaisir, rgPlaisir_sel)
# summary(rgPlaisir_sel)
# zPlaisir_sel = summary(rgPlaisir_sel)$coeff / summary(rgPlaisir_sel)$standard.errors
# pval_Plaisir_sel = 2 * (1 - pnorm(abs(zPlaisir_sel), 0, 1))
# pval_Plaisir_sel
# prPlaisir <- predict(rgPlaisir_sel)
# mcPlaisir <- table(PlaisirMod, prPlaisir)
# tauxPlaisir <- 1 - sum(diag(mcPlaisir))/sum(mcPlaisir)
# tauxPlaisir
# 
# # **************************************************
# # Régression pour Neophobie
# rDataNeophobie <- cbind(NeophobieMod, X_impute1)
# rDataNeophobie <- subset(rDataNeophobie, select = -c(Neophobie))
# 
# # mise en oeuvre de la regression
# library(nnet)
# rgNeophobie <- multinom(NeophobieMod~., data = rDataNeophobie)
# summary(rgNeophobie)
# zNeophobie = summary(rgNeophobie)$coeff / summary(rgNeophobie)$standard.errors
# pval_Neophobie = 2 * (1 - pnorm(abs(zNeophobie), 0, 1))
# pval_Neophobie
# 
# # Sélection de variables
# library(stats)
# rgNeophobie_sel <- step(rgNeophobie, k = 2, direction = 'both')
# # Comparaison de modèle avant et après sélection de variables
# anova(rgNeophobie, rgNeophobie_sel)
# summary(rgNeophobie_sel)
# zNeophobie_sel = summary(rgNeophobie_sel)$coeff / summary(
#   rgNeophobie_sel)$standard.errors
# pval_Neophobie_sel = 2 * (1 - pnorm(abs(zNeophobie_sel), 0, 1))
# pval_Neophobie_sel
# prNeophobie <- predict(rgNeophobie_sel)
# mcNeophobie <- table(NeophobieMod, prNeophobie)
# tauxNeophobie <- 1 - sum(diag(mcNeophobie))/sum(mcNeophobie)
# tauxNeophobie
# 
# # **************************************************
# # Régression pour Sélectivité
# rDataSelectivite <- cbind(SelectiviteMod, X_impute1)
# rDataSelectivite <- subset(rDataSelectivite, select = -c(Selectivite))
# 
# # mise en oeuvre de la regression
# library(nnet)
# rgSelectivite <- multinom(SelectiviteMod~., data = rDataSelectivite)
# summary(rgSelectivite)
# zSelectivite = summary(rgSelectivite)$coeff / summary(rgSelectivite)$standard.errors
# pval_Selectivite = 2 * (1 - pnorm(abs(zSelectivite), 0, 1))
# pval_Selectivite
# 
# # Sélection de variables
# library(stats)
# rgSelectivite_sel <- step(rgSelectivite, k = 2, direction = 'both')
# # Comparaison de modèle avant et après sélection de variables
# anova(rgSelectivite, rgSelectivite_sel)
# summary(rgSelectivite_sel)
# zSelectivite_sel = summary(rgSelectivite_sel)$coeff / summary(
#   rgSelectivite_sel)$standard.errors
# pval_Selectivite_sel = 2 * (1 - pnorm(abs(zSelectivite_sel), 0, 1))
# pval_Selectivite_sel
# prSelectivite <- predict(rgSelectivite_sel)
# mcSelectivite <- table(SelectiviteMod, prSelectivite)
# tauxSelectivite <- 1 - sum(diag(mcSelectivite))/sum(mcSelectivite)
# tauxSelectivite

