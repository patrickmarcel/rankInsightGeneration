library ( purrr )
library (moments)
# courbe roc 4 l i b r a r y (ROCR)
# pour arbre de d ́ecision 6 library ( rpart )
library ( rpart . plot )
# pour regression log
library(leaps)


N = 10000
n = c(10, 50, 100, 200, 500, 1000)
mu test = c(2, 5, 10, 100)
distri = c(”norm” , ”norm2” , ”unif” , ”khi” , ”exp”)
tests = expand.grid(distri , distri , n, n, mu test) 6 Etype I = rep(0, times = nrow(tests))
Skew1 = rep(0, times = nrow(tests))
Skew2 = rep(0, times = nrow(tests))
for(i in 1:nrow(tests)){ type 1=0
skew1 = 0
skew2 = 0
mu= tests[i, 5] for(k in 1:N){
e1 = 0 e2 = 0
# echant 1
i f ( t e s t s [ i , 1 ] == ” n o r m ” ) {
e1 = rnorm(tests[i, 3], mu, 1)
}
else if(tests[i,1]==”norm2”){
e1 = rnorm(tests[i, 3], mu, 20) }
else if(tests[i ,1] == ”unif”){
e1=runif(tests[i, 3],min=mu−1,max=mu+1)
}
else if(tests[i,1] == ”khi”){
e1 = rchisq(tests[i, 3], mu) }
else if(tests[i,1]==”exp”){ e1 = rexp(tests[i, 3], 1/mu)
}
# echant 2
if(tests[i, 2]==”norm”){
e2 = rnorm(tests[i, 4], mu, 1)
}
else if (tests[i, 2]==”norm2”){
e2 = rnorm(tests[i, 4], mu, 20) }
else if(tests[i, 2] == ”unif”){
e2=runif(tests[i, 4],min=mu−1,max=mu+1)
}
else if(tests[i, 2]==”khi”){
e2 = rchisq(tests[i, 4], mu) }
else if(tests[i, 2]==”exp”){ e2 = rexp(tests[i, 4], 1/mu)
}
type 1 = type 1 + as.numeric(t.test(e1, e2)$p.value<0.05) skew1 = skew1 + skewness(e1)
skew2 = skew2 + skewness(e2)
}
EtypeI[i]=type1/N
Skew1[i] = skew1/N
Skew2[i] = skew2/N }


data = data.frame(type I error = Etype I, n1 = tests[,3], n2 = tests[,4], skew1 =
Skew1 , skew2 = Skew2 )
# cr ́eation d’ indicateurs statistiques
# variable d’ interet
data$is out = as.numeric(data$type I error >0.075)
data$is out = as.factor(data$is out) 6 # statistiques
data$n1n2fois = data$n1∗data$n2
data$n1n2plus = data$n1 + data$n2
data$skew1skew2fois = data$skew1∗data$skew2
data$skew1skew2plus = data$skew1 + data$skew2
data$skew1skew2diff abs = abs(data$skew1 − data$skew2)
data$somme tot = (data$skew1/data$n1) + (data$skew2/data$n2)
data$diff tot abs = abs((data$skew1/data$n1) − (data$skew2/data$n2)) data$somme tot carre = data$skew1/(data$n1)ˆ2 + data$skew2/(data$n2)ˆ2


# r ́e ́equilibrage des 0−1 dans l ’ ́echantillon d’ apprentissage 
set . seed (1)
id eq = sample(which(data$is out==1), size = as.integer(0.8∗length(which(data$is out ==1)) ) )
id eq = c(id eq, sample(x = which(data$is out==0), size = length(id eq)))


# construction arbre
tree = rpart(is out  ̃ n1n2fois + n1n2plus + skew1skew2fois+ skew1skew2plus + skew1skew2diff abs + somme tot + diff tot abs + somme tot carre , data = data[id eq
,] , method = ”class”, cp = 10ˆ−10) 
# affichage de l ’ arbre construit
rpart . plot ( tree )
# fit sur apprentissage
# plot des donn ́ees en binaire
plot(data[id eq,]$diff tot abs, data[id eq,]$is out, ylab = ”hors limite (sup)”, xlab
= ”statistique”, main = ”Seuil sur  ́echantillon d’apprentissage”) 8 abline (v = 0.049 , col = ”red”)
# plot sur vraies donn ́ees
plot(data[id eq,]$diff tot abs, data[id eq,]$type I error , ylab = ”erreur type xlab = ”statistique”, main = ”Seuil sur  ́echantillon d’apprentissage”)
abline (v = 0.049 , col = ”red”)
abline (h = 0.075)
# m ́e triques apprentissage
pred tree = predict(tree , newdata = data[id eq,] , type = ”class”)
tp tree = sum(which(data[id eq,]$is out == 1)
fp tree = sum(which(data[id eq,]$is out == 0)
tn tree = sum(which(data[id eq,]$is out == 0)
fn tree = sum(which(data[id eq,]$is out == 1)
confusion matrix tree = matrix(c(tp tree , fp tree , fn tree , tn tree) , nrow = 2, byrow
= TRUE)
rownames(confusion matrix tree) = c(”predict 1”, ”predict 0”) colnames(confusion matrix tree) = c(”1”, ”0”)
print ( confusion matrix tree )
accuracy tree = (tp tree + tn tree)/length(pred tree) precision tree = tp tree/(tp tree + fp tree)
recall tree = tp tree/(tp tree + fn tree)
f2 score tree = 5∗(recall tree∗precision tree)/(recall tree + 4∗precision tree) cat(”\n accuracy: ”)
cat(accuracy tree)
cat(”\n precision: ”)
cat(precision tree)
cat(”\n recall : ”)
cat(recall tree)
cat(”\n f2 score: ”)
cat(f2 score tree)
# pr ́e dictions sur test
# plot des donn ́ees en binaire
plot(data[−id eq,]$diff tot abs, data[−id eq,]$is out, ylab = ”hors limite (sup)”, xlab = ”statistique”, main = ”Seuil sur  ́echantillon test”)
abline (v = 0.049 , col = ”red”) # plot sur vraies donn ́ees
plot(data[−id eq,]$diff tot abs, data[−id eq,]$type I error , ylab = ”erreur type 1”, xlab = ”statistique”, main = ”Seuil sur  ́echantillon test”)
abline (v = 0.049 , col = ”red”)
abline (h = 0.075) #m ́etriques test
data test = data[−id eq,]
pred tree = predict(tree , newdata = data test , type = ”class”) tp tree = sum(which(data[−id eq,]$is out == 1) %in% which(pred fp tree = sum(which(data[−id eq ,]$is out == 0) %in% which(pred tn tree = sum(which(data[−id eq,]$is out == 0) %in% which(pred fn tree = sum(which(data[−id eq ,]$is out == 1) %in% which(pred confusion matrix tree = matrix(c(tp tree , fp tree , fn tree , tn
= TRUE)
rownames(confusion matrix tree) = c(”predict 1”, ”predict 0”) colnames(confusion matrix tree) = c(”1”, ”0”)
print ( confusion matrix tree )
accuracy tree = (tp tree + tn tree)/length(pred tree) precision tree = tp tree/(tp tree + fp tree)
tree==1))
tree==1))
tree==0))
tree==0)) tree),nrow=2,byrow
recall tree = tp tree/(tp tree + fn tree)
f2 score tree = 5∗(recall tree∗precision tree)/(recall tree + 4∗precision tree) cat(”\n accuracy: ”)
cat(accuracy tree)
cat(”\n precision: ”)
cat(precision tree)
cat(”\n recall : ”)
cat(recall tree)
cat(”\n f2 score: ”)
cat(f2 score tree)


# s ́election backward
mod. sel back = regsubsets ( is out  ̃ 0 + n1n2fois + n1n2plus + skew1skew2plus +
skew1skew2diff abs + somme tot + diff tot abs + somme tot carre , data = data[id eq
,] , method = ”backward”)
summary(mod. sel back)
# construction mod`ele
mod log = glm(is out  ̃ 0 + diff tot abs, data[id eq,], family = binomial(link = logit)
)
summary(mod log)
# choix du seuil de proba
pred = predict(mod log, data[id eq,], type = ”link”)
pred = prediction(pred, data$is out[id eq]) roc = performance(pred,”tnr”,”fnr”) plot(roc, colorize=T, lwd=2)
abline(a = 0, b = 1)
threshold = roc@alpha.values[[1]][which(roc@x.values[[1]]<=0)[1]]# on limite le taux de faux n ́egatifs `a 0%
print(threshold)
# seuil sur statistique directement
threshold stat = threshold/mod log$coefficients [[1]]
threshold stat
# fit sur apprentissage
# plot des donn ́ees en binaire
plot(predict(mod log, newdata = data[id eq,], type = ”link”), data[id eq,]$is out, ylab = ”hors limite (sup)”, xlab = ”statistique”, main = ”Seuil sur  ́echantillon d’ apprentissage”)
abline (v = threshold , col = ”red”) # plot sur vraies donn ́ees
plot(predict(mod log, newdata = data[id eq,], type = ”link”), data[id eq,]$type I
error , ylab = ”erreur type 1”, xlab = ”statistique”, main = ”Seuil sur  ́echantillon d’ apprentissage”)
abline (v = threshold , col = ”red”) 
abline (h = 0.075)

# m ́e triques apprentissage
pred log = predict(mod log, newdata = data[id eq,], type = ”link”)>=threshold tp log = sum(which(data[id eq,]$is out == 1) %in% which(pred log == 1))
fp log = sum(which(data[id eq,]$is out == 0) %in% which(pred log == 1))
tn log = sum(which(data[id eq,]$is out == 0) %in% which(pred log == 0))
fn log = sum(which(data[id eq,]$is out == 1) %in% which(pred log == 0))
confusion matrix log = matrix(c(tp log , fp log , fn log , tn log) , nrow = 2, byrow =
TRUE)
rownames(confusion matrix log) = c(”predict 1”, ”predict 0”) colnames(confusion matrix log) = c(”1”, ”0”)
print ( confusion matrix log )
accuracy log = (tp log + tn log)/length(pred log)
precision log = tp log/(tp log + fp log)
recall log = tp log/(tp log + fn log)
f2 score log = 5∗(recall log∗precision log)/(recall log + 4∗precision log) cat(”\n accuracy: ”)
cat ( accuracy log )
cat(”\n precision: ”)
cat(precision log)
cat(”\n recall : ”)
cat(recall log)
cat(”\n f2 score: ”)
cat(f2 score log)
# pr ́e dictions sur test
# plot des donn ́ees en binaire
plot(predict(mod log, newdata = data[−id eq,], type = ”link”), data[−id eq,]$is out,
ylab = ”hors limite (sup)”, xlab = ”statistique”, main = ”Seuil sur  ́echantillon test”)
abline (v = threshold , col = ”red”) # plot sur vraies donn ́ees
plot(predict(mod log, newdata = data[−id eq,], type = ”link”), data[−id eq,]$type I error , ylab = ”erreur type 1”, xlab = ”statistique”, main = ”Seuil sur  ́echantillon test”)
abline (v = threshold , col = ”red”)
abline (h = 0.075) #m ́etriques test
data test = data[−id eq,]
pred log = predict(mod log , newdata = data test , type = ”link”)>=threshold
tp log = sum(which(data[−id eq,]$is out == 1) %in% which(pred log == 1))
fp log = sum(which(data[−id eq,]$is out == 0) %in% which(pred log == 1))
tn log = sum(which(data[−id eq,]$is out == 0) %in% which(pred log == 0))
fn log = sum(which(data[−id eq,]$is out == 1) %in% which(pred log == 0))
confusion matrix log = matrix(c(tp log , fp log , fn log , tn log) , nrow = 2, byrow =
TRUE)
rownames(confusion matrix log) = c(”predict 1”, ”predict 0”) colnames(confusion matrix log) = c(”1”, ”0”)
print ( confusion matrix log )
accuracy log = (tp log + tn log)/length(pred log)
precision log = tp log/(tp log + fp log)
recall log = tp log/(tp log + fn log)
f2 score log = 5∗(recall log∗precision log)/(recall log + 4∗precision log) cat(”\n accuracy: ”)
cat ( accuracy log )
cat(”\n precision: ”)
cat(precision log)
cat(”\n recall : ”)
cat(recall log)
cat(”\n f2 score: ”)
cat(f2 score log)