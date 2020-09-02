# Sandra Tretter
# 14.08.2020

###############################
###   FISHER'S EXACT TEST   ###
###############################

#	content:	0.	required packages

#		      	1.	Pogonomyrmex worker against Camponotus male

#           2. Paired, P. califonicus

#           3. Paraffin vs. non-paraffin

#-----------------------------------------------------------------------------------------------------------------------------#

####	0.	required packages ####

#install.packages("Exact")		# remember: small and capital letters matter!!!
library(Exact)				# if you want to use the Barnard test

#install.packages("exact2x2")		# for an alternative approach, apparently the confidence intervalls are better
library("exact2x2")			# calculated than in the built in fisher-exact-test


#-----------------------------------------------------------------------------------------------------------------------------#

####	1. 	Pogonomyrmex worker against Camponotus male ####

#               gestochen (+)   nicht gestochen (=)   Summe der Versuche
# Pogo           2                   62                      64
# Camponotus     8                    4                      12


# make up the matrix (and call it "Camponotus"):


Camponotus <-			matrix(c(64, 12, 2, 8),					# make a matrix with these values (rows first)
                       nrow = 2,							# and "number of rows" (nr) =2
                       dimnames = list(c("Pogonomyrmex", "Camponotus"),	# give names to the rows (first), 
                                       "Agression"= c("Replicates", "Stinging/C-Posture")))	# then columns (and give a title for the column and rows, if you want)


Camponotus						# have a look at the data

fisher.test(Camponotus)				# run the Fisher exact test 
# result:(p-value = 7.949e-05)

chisq.test(Camponotus)					# run the chi square test
# result: X-squared = 16.976, df = 1, p-value = 3.786e-05

exact.test(Camponotus)       # run the Bernard test
# p-value = 2.951e-05

#-----------------------------------------------------------------------------------------------------------------------------#

#### 2. Paired, P. californicus ####

## own vs. foreign colonies

#                Own     Foreign   
# Observed       4       21                                 
# Expected       12      13                                  

Biting1 <- matrix (c(4,12,21,13), nrow = 2, dimnames = list(c("Observed","Expected"), c("Own","Foreign"))) # 20 bites haplo + 5x first biter in HP
Biting1

fisher.test(Biting1)                  
# result:(p-value = 0.03216)

Biting2 <- matrix (c(4,12.5,21,12.5), nrow = 2, dimnames = list(c("Observed","Expected"), c("Own","Foreign"))) # 20 bites haplo + 5x first biter in HP
Biting2

chisq.test(Biting2)
# result: X-squared = 5.0882, df = 1, p-value = 0.02409

exact.test(Biting1)
# p-value = 0.01637

## haplo vs. pleo

#                Haplo   Pleo   
# Observed       25       8                                 
# Expected       17      16 

Biting3 <- matrix (c(25, 17, 8, 16), nrow = 2, dimnames = list(c("Observed","Expected"), c("Haplo","Pleo"))) # 20 bites haplo + 5x first biter in HP
Biting3

fisher.test(Biting3)                  
# result:(p-value = 0.07224)

exact.test(Biting3)
# p-value = 0.04591

Biting4 <- matrix (c(25, 16.5, 8, 16.5), nrow = 2, dimnames = list(c("Observed","Expected"), c("Haplo","Pleo"))) # 20 bites haplo + 5x first biter in HP
Biting4

chisq.test(Biting4)
# result: X-squared = 3.6513, df = 1, p-value = 0.05602


# Binomial Test

biting5 <- c(25,8)
binom.test(biting5)
# result: p-value = 0.004551


#------------------------------------------------------------------------------------------------------------------------------#

#### 3. Paraffin / non-Paraffin ####

#                with Paraffin   without Paraffin   
# Observed       253             461                                 
# Expected       357             367 

Paraffin <- matrix(c(253,357,461,357), nrow = 2, dimnames = list(c("Observed","Expected"), c("with Paraffin","without Paraffin"))) # count of interactions
Paraffin

fisher.test(Paraffin)
# result: p-value = 3.381e-08

chisq.test(Paraffin)
# result: X-squared = 30.361, df = 1, p-value = 3.586e-08

exact.test(Paraffin)
# p-value = 2.647e-08