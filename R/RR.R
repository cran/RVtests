RR <-
function(x, y, scale = FALSE, lambda = 1, npermutation = 100, npermutation.max, min.nonsignificant.counts)
{
#RR: ridge regression
#x=UDV': svd(x)
#always centralizing x and y
x<- scale(x, scale = scale)
y<- y-mean(y)

if (missing(min.nonsignificant.counts)) min.nonsignificant.counts<- 10  #?10, ?20
if (missing(npermutation.max)) npermutation.max<- npermutation

UDV<- svd(x)
	d<- UDV$d
	u<- UDV$u
	v<- UDV$v
	indx<- (abs(d) > 1e-8)  #remove zero eigenvalue.
	d<- d[indx]
	u<- u[,indx, drop=FALSE]
	v<- v[,indx, drop=FALSE]

#test score and pvalue (nominal pvalue)

yp<- y
uy<- t(u)%*%yp
sco<- rep(NA, length(lambda))
names(sco)<- paste("RR", lambda,sep="")
for (i in 1:length(lambda)){
	lam<- lambda[i]
	ye<- u%*%((d^2/(d^2+lam))*uy)
	sco[i]<- cor(yp, ye)
	}
testscore<- abs(sco)
testpvalue<- NULL
#testscore<- sco*sqrt( (length(yp)-2)/(1-sco^2) ) #t-distribution with df of n-2. not correct??? cor(y, Ay)?
#testpvalue<- 2*pt(abs(testscore), df = length(yp)-2, lower.tail = FALSE)

#permutation pvalue (empirical pvalue)
#permutation
permpvalue<- NULL
counts<- NULL
jth<- 0

if (npermutation >= 1) {
counts<- rep(0, length(lambda))
while ((jth < npermutation) | ((jth < npermutation.max) & (min(counts) < min.nonsignificant.counts)))
{
jth<- jth + 1
	yp <- y[sample(1:length(y), replace=FALSE)]
	uy<- t(u)%*%yp
	sco<- rep(NA, length(lambda))
	#names(sco)<- paste("RR", lambda,sep="")
	for (i in 1:length(lambda)){
		lam<- lambda[i]
		ye<- u%*%((d^2/(d^2+lam))*uy)
		sco[i]<- cor(yp, ye)
		}
	permscore<- abs(sco)
	#permscore<- sco*sqrt( (length(yp)-2)/(1-sco^2) ) #t-distribution with df of n-2.???
	counts<- counts + (permscore >= testscore)
	}#while
	permpvalue<- (1+ counts)/(1 + jth)
}#if

list(score = testscore, nonsignificant.counts = counts, pvalue.empirical = permpvalue, pvalue.nominal = testpvalue, total.permutation = jth)
}#RR function

