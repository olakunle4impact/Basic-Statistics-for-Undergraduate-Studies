getwd()
Result<- data.frame(names=c("Joshua","Kezia","John","Sarah"),
          Department=c("Sales","Marketing","Admin","Research"),
          Salary=c(105000, 200050, 305000, 250000))
Result
variable1=c(1,2,3,4)
cat("variable1 is", "ade", variable1)
Even=c(2,4,6,8)
cat("The second variable", Even, "is a list of four consecutive even numbers")
print(ls(all.names=T))
print(ls(all.names = T))
2^2/3/((1^1/5)-(1/3))
2^2/3
(1/3)+((5/7)*((9/10)-1+(3/4)))
2^(2+2)-4+16^((-2)^(2.25-(1/4)))
log(243,base = 3)
exp(x=3)
4^3*(5-3)/(51-3+7)
5^2+3*50/8-3
seq(from=5,to=10,by=.5)
rep(x=2, times = 4)
rep(score=c(2,4,5), times=2)
rep(score=c(2,3,5), times=2)
rep(score=c(3,62,8.3), times=3)
new=c(2,3,4)
rep(new,3)

scores<-c(12,14,18,11,16)
scores[1]

kay <- list(c(T,F,F,T),"Olakunle",matrix(data = c(2:5), nrow = 2, ncol = 2))
kay
kay[[3]][2,2]
names(kay)<-c("Logical","String","Matrix")
kay$Matrix[2,1]
kay$String
kay$Numeric<-c(10:13)
names(kay)

a=c(1,2,3)
b=c(2,3,4)
rbind(a,b)

cbind(a,b)

matrix<-matrix(c(1:9), 3, 3)
matrix

matrix[3,2]
matrix[,2]
matrix[2:3,]

matrix[c(3,1,2),2:3]
matrix[-1,2:3]
matrix[-1,-2]
matrix[-1,-c(2,3)]
matrix[2,]<- c(3,4,5)
matrix
matrix[c(1,3),2]<-c(8,2)
matrix[c(1,3),2]
diag(x=matrix)
diag(x=matrix)<-rep(x=0, times = 3)
matrix
t(matrix)

diag(x=3)
a<-array(c(2,3,1),dim=c(2,3,3,3))
a
seq(2.5,10,0.5)
Names <- c("James", "Joshua", "kezia", "Kole", "John")
Department <- c("Sales", "Finance", "Marketing", "Sales", "Admin")
Income <- c(150000, 230000, 100050, 150000, 205000)
Employee <- data.frame(Names, Department, Income)
Employee


employee <- data.frame(names = c("James","Joshua","Kezia","Kole","John"),
                       Dept = c("Sales","Finance", "Marketing","Sales","Admin"),
                       income = c(150000, 230000, 100050, 150000, 205000))
employee

employee[3:5,2]
employee$names[3]
dim(employee)
newEmployee<-data.frame(names = "Joseph", 
                        Dept = factor("Admin",levels = levels(employee$Dept)), 
                        income = 350000)

employee<-rbind(employee,newEmployee)
employee
newEmployee

Gender <- c("Male","Male","Female","Male","Male","Male")
Gender <-factor(x=Gender, levels = c("Male","Female"))
Gender
employee <- cbind(employee,Gender)
employee

employee[employee$Gender=="Male", c("names","income")]
employee[employee$income>200000&employee$Gender=="Male", ]
typeof(employee)
class(employee)

Months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
Months <-factor(x=Months, levels = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"), ordered = T)
Months

ls("package:reshape")
search()

ls("package:dplyr")
environment(getwd)
library(caret)
library("dplyr")

cummean(c(1,2,2,2,3,4,5))
5/3

row.bind()
class.scores

seq(by = 0.2, to = 6, from = 3)
seq(b = 0.2, t = 6, f = 3)
seq(3,6,0.2)
seq(from = 3, 6, b = 0.2)
seq(to = 6, f = 3, 0.2)

args(c)

args(data.frame)
args(plot)

b <- 5
if(b<=5){
  b^2
}

vtr <- c(2,4,6,5,7,1,8)
mat <- matrix(vtr,3,3)
if(any(vtr-1)>8||mat[2,2]<=8) {
  cat("if true, condition is satisifed.. \n")
  mylist = list(aa = vtr, bb = mat)
  cat("a list with", length(mylist), "is generated." )
}

mylist$aa

seq(1,9,2)
vtr[seq(1,9,2)]=NA
vtr

mat[2,2]

if(any(vtr-1)>8){
  cat("This condition is true. \n")
  }else{
    cat("This condition is false")
  }

x<-4
y<-c(-4:4)
y==0
score<-ifelse(test = y==0, yes = NA, no = x/y)
score

for (variable in vector) {
  #put your code here
}
 
boat <- 2
for (variable in 3:7) {
   boat = boat +1
   cat("As", boat, "increases, so does,", variable,"\n")
}

var<-c(2:8)
for (i in var) {
  print(2^i)
}


kay <- 4
while(kay<3){
  kay <- kay + 1
}

scr <-function(){
  scr1<- 1
  scr2<- 1
  cat(scr1, ",", scr2, ",", sep = "")
  repeat{
    result<-scr1+scr2
    scr1<-scr2
    scr2<-result
    cat(scr2,",", sep = "")
    if(scr2>100){
      cat("Break right away")
      break
    }
  }
}

scr <-function(limit){
  scr1<- 1
  scr2<- 1
  cat(scr1, ",", scr2, ",", sep = "")
  repeat{
    result<-scr1+scr2
    scr1<-scr2
    scr2<-result
    cat(scr2,",", sep = "")
    if(scr2>limit){
      cat("Break right away")
      break
    }
  }
}

scr(300)



ls()

f <- function(x,y){
  x^2 + y /z
}
environment(f)
search()

f(2,3)
z=2

make.power <- function(n){
  pow <- function(x){
    x^n
  }
  pow
}

x=3
make.power(4)
pow(3)
cube<-make.power(3)
cube
square<-make.power(2)
square
cube(3)
square(7)
pentic<-make.power(4)
pentic
scr()
scr(limit = )

r<-function(x,y){
  x^3+ y /z
}
r(2,3)
z=2
3/2

y<- 8
k <- function(x){
  y<-3
  y^2 + g(x)
}

g<-function(x){
  x*y
}

k(3)
g(3)


a=3
k <- function(x){
  y<-4
  y^2 + g(z)
}
k(3)

ls(environment(k))
get("y",environment(k))


k <- function(x){
  y<-2
  y^2 + g(x)
}

y<- 8
k <- function(x){
  y<-3
  y^2
}
k(3)

g<-function(x){
  x*y
}

g(3)

bay <- function(a,b){
  a^3
}
bay(3)


?trees
print(trees)

set.seed(11223)
sample(trees$Height)

set.seed(12134)
sample(trees$Height, replace = F)

sample(c("Johnson","George","Kezia","Joshua","Samuel"),1)
a=5

set.seed(22312)
sample(trees$Height, size = 5, replace = T)

set.seed(13131)
sample(trees$Height, size = 5, replace = F)

set.seed(12441)
sample(c("Good","Bad"), size = 5, replace = F, prob = c(.75,.25))


trees[sample(nrow(trees),5), ]

library(dplyr)
sample_n(trees,5)

?split
?cut
?cut
cut(trees$Height, breaks = 6)

sample(trees$Height, size = 6, replace = F)
sample(trees$Height, size = 6, replace = F)
sample(trees$Height, size = 6, replace = F)
sample(trees$Height, size = 6, replace = F)
sample(trees$Height, size = 6, replace = F)
sample(trees$Height, size = 6, replace = F)
sample(trees$Height, size = 6, replace = F)


set.seed(14511)
start=sample(1:5,1)
start

chickwts

split(chickwts$weight,chickwts$feed)
?split

chickwts$weight[chickwts$feed=="sunflower"]
table(chickwts$feed)
attach(chickwts)

feedC=subset(chickwts,feed=="casein")       
feedC


feedC[sample(nrow(chickwts$feed), 3), ]

mathScore<-c(12,17,14,11,19,14,16,11,10,18)
length(mathScore)
mean(mathScore)

example<-c(68,70,72,41,50)

geoMean <- function(dataset){
  prod(dataset)^(1/length(dataset))
}

geoMean(example)
geoMean(mathScore)

chemScore<-c(13,12,11,15,17,19,20,15,16,100)
MeanChem<-mean(chemScore)

harMean<-function(dataset){
  (1/length(dataset)*sum(1/dataset))^(-1)
}

harMean(mathScore)
mathScore
median(mathScore)
?match

mathScore<-sort(mathScore)
median(mathScore)

mathQuantile<-quantile(mathScore, probs = c(seq(0,1,0.20)))
mathQuantile

summary(mathScore)

tagMode<-function(dataset)
{
  UniDataSet <- unique(dataset)
  UniDataSet[which.max(tabulate(match(dataset,UniDataSet)))]
}

mode(mathScore)
?which.max
unique(mathScore)
?tabulate
tabulate(chemScore)

tagMode(mathScore)

?summary
range(mathScore)
mathRange<-max(mathScore)-min(mathScore)
mathRange

IQR(mathScore)

str(mathScore)
mad(mathScore)
cumsum(mathScore)
cummax(mathScore)


library(help="stats")

cube <- function(x, n) {
  x^3
}
cube(3)
x <- 1:10
if(x > 5) {
  x <- 0
}

f <- function(x) {
  g <- function(y) {
    y + z
  }
  z <- 4
  x + g(x)
}

z <- 10
f(3)
x <- 5
y <- if(x < 3) {
  NA
} else {
  10
}

h <- function(x, y = NULL, d = 3L) {
  z <- cbind(x, d)
  if(!is.null(y))
    z <- z + y
  else
    z <- z + f
  g <- x + y / z
  if(d == 3L)
    return(g)
  g <- g + 10
  g
}

pollutantmean<-function(directory, pollutant, id = 1:332){
  ## "directory is a character vector of length 1 indicating
  ## the location of the CSV files
}

complete<-function(directory, id = 1:332){}
cor<-function(directory, threshold = 0){}


source("complete.R")
complete("specdata", 1)
##   id nobs
## 1  1  117
complete("specdata", c(2, 4, 8, 10, 12))
##   id nobs
## 1  2 1041
## 2  4  474
## 3  8  192
## 4 10  148
## 5 12   96
complete("specdata", 30:25)
##   id nobs
## 1 30  932

source("pollutantmean.R")
pollutantmean("specdata", "sulfate", 1:10)
## [1] 4.064128
pollutantmean("specdata", "nitrate", 70:72)
## [1] 1.706047
pollutantmean("specdata", "nitrate", 23)
## [1] 1.280833


## 2 29  711
## 3 28  475
## 4 27  338
## 5 26  586
## 6 25  463
complete("specdata", 3)
##   id nobs
## 1  3  243


source("corr.R")
source("complete.R")
cr <- corr("specdata", 150)
head(cr)
## [1] -0.01895754 -0.14051254 -0.04389737 -0.06815956 -0.12350667 -0.07588814
summary(cr)
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## -0.21057 -0.04999  0.09463  0.12525  0.26844  0.76313
cr <- corr("specdata", 400)
head(cr)
## [1] -0.01895754 -0.04389737 -0.06815956 -0.07588814  0.76312884 -0.15782860
summary(cr)
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## -0.17623 -0.03109  0.10021  0.13969  0.26849  0.76313
cr <- corr("specdata", 5000)
summary(cr)
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 
length(cr)
## [1] 0
cr <- corr("specdata")
summary(cr)
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## -1.00000 -0.05282  0.10718  0.13684  0.27831  1.00000
length(cr)
## [1] 323

sd(mathScore)
library(help="stats")
environment(var)
var(mathScore)

sqrt(var(mathScore))

SEE<- sd(mathScore)/sqrt(length(mathScore))
print(SEE)


set.seed(12334)
normal<-rnorm(20,0,1)
mean(normal)
sd(normal)

library(moments)
skewness(normal)

kurtosis(normal)

data <- c(0.7,-0.6,0.1,-0.1,-1.0,0.4,0.3,-1.7,1.1,4.0)
plot(data,rep(0,10),axes=T,ylab="",cex=2,cex.axis=1.5,cex.lab=1.5)
abline(h=0,col="gray",lty=2)
arrows(3,0.3,3.9,0.1,lwd=2)
text(3,0.45,labels="Outlier?",cex=2)

#To get bivariate example;
bee <- c(0.1,0.3,1.3,0.6,0.2,-1.7,0.8,0.9,-0.8,-1.0)
boo <- c(-0.3,0.9,2.8,2.3,1.2,-4.1,-0.4,4.1,-2.3,-100.0)
plot(bar,baz,axes=T,cex=2,cex.axis=1,cex.lab=1.5)
arrows(-0.5,-80,-0.94,-97,lwd=2)
text(-0.45,-74,labels="Outliers?",cex=3)

library(stats)
ls("package:stats")

colSums(trees)
colMeans(trees)
rowSums(trees)
rowMeans(trees)


cumsum(mathScore)
cummax(mathScore)
cummin(mathScore)
cumprod(mathScore)

department = c("EMT", "WMA", "AQFM", "FWM")
department[4]

c("Joshua", "James", "John", "Jordan") -> Studentnames
Studentnames[2]

scores<-c(2,5,7,3,NaN,2,6,4,NaN,1)
is.na(scores)
table(is.na(scores))

!is.na(scores)
table(!is.na(scores))

mathScore
?cummax
cummax(v)
v=c(1,3,5,6,7,2,8,5,6,9)

mtree<-c(25.46,31.51,17.04,30.10,38.40,31.19,37.56,41.70,21.01,41.83)
n<-length(mtree)
print(n)
mtree.mean<-mean(mtree)
print(mtree.mean)
mtree.sd<-sd(mtree)
print(mtree.sd)
mtree.se<- mtree.sd/sqrt(n)
print(mtree.se)
mtree.test<- (mtree.mean - 28)/mtree.se
pt(mtree.test, df=n-1)

m.test<-t.test(x=mtree,mu=28,alternative="less")
m.test

?pt
mtree.mean

consumption<-read.csv("Standard dev paper.csv")
consumption.means<-tapply(consumption$Energy,INDEX = consumption$Type,FUN=mean)
boxplot(consumption$Energy~consumption$Type)
points(1:4,consumption.means,pch=4,cex=1.5)
consumption.sd<-tapply(consumption$Energy,INDEX = consumption$Type,FUN=sd)
max(consumption.sd)/min(consumption.sd)
consumption.meancen<-consumption$Energy-consumption.means[as.numeric(consumption$Type)]
qqnorm(consumption.meancen,main="Normal QQ plot of residuals")
qqline(consumption.meancen)
attach(consumption)
consumption.anova<-aov(Energy~Type,data=consumption)
summary(consumption.anova)

kruskal.test(Energy~Type,data = consumption)

data<-read.csv("practice.csv")
data
names(data)
str(data)
data.result<-aov(yield~cultivars,data=data)
data.result
summary(data.result)
options(show.signif.stars=F)

plot(data.result$fitted.values,data.result$residuals,
     main = "Residuals Vs Fitted", pch =20)
abline(h=0, lty=2)


n1<-sum(data.result$model$cultivars =="A")
n2<-sum(data.result$model$cultivars =="B")
n3<-sum(data.result$model$cultivars =="C")

sum(data$cultivars=="A")
sum(data$cultivars=="B")
sum(data$cultivars=="C")
as.integer(sum(data$cultivars=="A"))
MSE<-sqrt(sum((data.result$residuals)^2)/data.result$df.residual)
tcrit <- qt(0.025, data.result$df.residual, lower.tail = F)
LSD<- tcrit*sqrt((2*MSE)/n1)
print(LSD)

ab= mean(data$yield[data$cultivars=="A"])
bm= mean(data$yield[data$cultivars=="B"])
cm= mean(data$yield[data$cultivars=="C"])
print(ab)
print(bm)
print(cm)

data.tukey<-TukeyHSD(data.result, ordered = T)
print(data.tukey)
plot(data.tukey)

practice1<-read.csv("practice1.csv")
practice1
## Data manipulation. log transformation ##
transform<-log10(practice1$mp) ## you can also use "log(practice1$mp,10) i.e base 10 ##
## log(practice$mp) will produce the natural logarithm ##
transform
log(practice1$mp,10)
dataP<-data.frame(nlevel=practice1$nlevel,mp=transform) ## renaming transform ##
dataP<-data.frame(nlevel=c(practice1$nlevel),mp=round((c(transform)),2))
dataP<-data.frame(nlevel=c(as.character(practice1$nlevel)),mp=round((c(transform)),2))
dataP
names(dataP)
dataP$mp[dataP$nlevel=="A"] ## mp for nlevel A only ##
output<-aov(mp~nlevel,data=dataP)
anova
summary(output)


math_score<-c(18,12,15,11,15,17,17,16,18,18)
barplot(math_score, main = "Students Scores in Mathematic", xlab = "scores",
        ylab = "Department", col = "darkblue", density= 10,
        border = "red", horiz = F)
dev.new()
barplot(math_score, main = "Students Scores in Mathematic", xlab = "scores",
        ylab = "Department", col = "darkblue", density= 10,
        border = "red", horiz = T)
dev.new(2)
hist(math_score)

colerm<-matrix(c(15,11,18,14,12,15,18,12),2,4,
              dimnames=list(c("EMT","FWM"),c("P_score","M_score","C_score","B_score")))
barplot(colerm, main = "First Year Score of Colerm Student", xlab = "Subject",
        ylab = "Score", col = c("darkblue","darkgreen"))
        legend("topleft", c("EMT","FWM"), fill = c("darkblue","darkgreen"))
print(colerm)
?legend
hist(math_score, main = "Hisotgram of Math_score", xlab = "Math_score", xlim = c(10,20),
     col = "brown")

hist<-hist(math_score,main = "Hisotgram of Math_score", xlab = "Math_score",
           col="red")
text(hist$mids,hist$counts,labels=hist$counts, adj=c(0.5,-0.5))
print(hist)

library(stats)

ls("package:MASS")
polygon(math_score)

scores<-c(18,10,112,13,13,14,16,18,17,18)
library(MASS)
frequency.polygon(scores)


getwd()

latin<-read.csv("latin.csv")
names(latin)
print(latin)
class(latin)
str(latin)
latin$Row<-factor(latin$Row)
latin$Col<-factor(latin$Col)
latin$Drug<-factor(latin$Drug)
str(latin)

mod3<-aov(latin$SodiumLevel~latin$Row+latin$Col+latin$Drug, data = latin)
summary(mod3)
anova(mod3)
names(mod3)
plot(mod3$fitted.values,mod3$residuals, main = "Residuals Vs Fitted", pch=20)
abline(h=0, lty= 2)

plot(mod3$model$`latin$Drug`,mod3$residuals)

model.tables(mod3,"effects")
replications(latin$SodiumLevel~latin$Row+latin$Col+latin$Drug, data = latin)
!is.list(replications(latin$SodiumLevel~latin$Row+latin$Col+latin$Drug, data = latin))
?cor

?trees
library(help="datasets")
names(trees)
str(trees)

cor(trees$Girth,trees$Volume, method = "pearson" , use = "complete.obs")
cor(trees$Girth,trees$Volume, method = "kendall" , use = "complete.obs")
cor(trees$Girth,trees$Volume, method = "spearman" , use = "complete.obs")

cor(trees$Girth,trees$Height)
cor(trees$Height,trees$Volume)
cor(trees$Volume,trees$Height)

plot(trees$Volume,trees$Girth, xlab = "Height", ylab = "Volume")
slr<-lm(trees$Volume~trees$Height,data = trees)
print(slr)
abline(slr1,lwd=1)
slr1<-lm(trees$Volume~trees$Girth,data = trees)
names(slr)

summary(slr)
names(summary(slr))
summary(slr)$sigma
trees
slr
detach(trees)
data<-data.frame(Height=c(35:40))
data
predV<-predict(slr,newdata = data, interval = "confidence", level = 0.95)
predV
library(stats)
survey
length(survey$Wr.Hnd)
ls("package:stats")
library(MASS)
survey
trees

?predict

attach(trees)
model<-lm(Volume~Height, data = trees)
predict(model, newdata = data.frame(Height=76), interval = "confidence",
        level = 0.95)

predict(model, newdata = data, interval = "confidence",
        level = 0.95)

plot(Volume~Height)
trees

bal<-data.frame(Wr.Hnd=c(14.5,24))




sol<-lm(Height~Wr.Hnd, data = survey)
print(sol)
mol<-predict(sol,newdata = bal,interval = "confidence", level = 0.95)
mol

class(survey)
class(trees)
trees

attach(trees)
mlr<-lm(Volume~Girth+Height, data = trees)
summary(mlr)

predict(mlr, newdata = data, interval = "prediction", level = 0.99)

data<-data.frame(Girth =c(11:15), Height =c(76:80))

mlrP<-lm(Volume~Girth+I(Girth^2), data = trees)
summary(mlrP)

mlrT<- lm(Volume~Girth+Height+Girth:Height, data = trees)
summary(mlrT)

mlrC<- lm(Volume~Girth*Height, data = trees)
summary(mlrC)


nlr<-read.csv("nonlinear.csv")
nlr
attach(nlr)
plot(X,Y, main = "Scatterplot showing relationship between X and Y")
model<-lm(log(Y)~log(X),data=nlr)
summary(model)
ls("package:lattice")
library(ggplo)
a=table(trees$Height)
frequency.polygon(a)
scores<-c(10,17,19,11,14,16,18,15,18,13)
simple.freqpoly(scores)

hist(scores)
frequency.polygon(table(scores))
table(scores)

plot(scores, type = "l", main = "Line graph of Scores")

boxplot(trees, main = "Boxplot for trees", horizontal = T, notch = T,
        col = c("red","green","blue"))

stripchart(obj, main = "Stripchart for Height and Volume", method = "jitter",
           col = c("red","green"))
Height<-trees$Height
Girth<-trees$Girth
Volume<-trees$Volume
obj<-list("Height"=Height, "Volume"=Volume)


stripchart(trees$Volume~trees$Girth, main ="Stripchart for Volume~Girth", col = "red")
attach(trees)
plot(Volume~Girth, main = "scatterplot of Volume~Girth")
plot(Girth,Volume)


scores<-c(16,20,14,18)
dept<-c("EMT","FWM","AQFM","WMA")
pie(scores,labels = dept, main = "Pie Chart of Colerm Students Scores")
percent<-round(scores/sum(scores)*100)
dept<-paste(dept,percent)
dept<-paste(dept,"%", sep = "")
pie(scores, labels = dept, col = rainbow(length(dept)),
    main = "Pie Chart of Colerm Students Scores")
