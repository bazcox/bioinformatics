# T = 2; C = 3 conversion

b0 <-	0.813
b1 <-	-7.104E-5
b2 <-	2.845E-4
b3 <-	-1.992E-4
b4 <-	-2.062E-4
b5 <-	-1.089E-5
b6 <-	3.588E-4
b7 <-	2.978E-5
b8 <-	2.477E-4
b9 <-	-1.737E-4
b10 <-	-2.277E-4
b11 <-	-6.214E-5
b12 <-	-1.503E-4
b13 <-	-9.380E-4
b14 <-	2.230E-4
b15 <-	-1.306E-4
b16 <-	-0.531
b17 <-	0.239
b18 <-	-0.003
b19 <-	0.025
b20 <-	-0.01
b21 <-	-5.631E-4
b22 <-	-0.001
b23 <-	-1.940E-4
b24 <-	-1.715E-4
b25 <-	4.414E-4
b26 <-	-0.044
b27 <-	-0.084
b28 <-	-0.007

# chnage filename accordingly
JASPfile <- read.csv("GWAS_75testset_all_afterBH_probe5_v2_converted_dummycode.csv")
JASPfile <- na.omit(JASPfile)

# change variable names (rsxxxx_cat) & bxx accordingly
# predictions <- (b0 + b1*JASPfile$rs1034420_cat + b2*JASPfile$rs78578442_cat
#                 + b3*JASPfile$rs2003734_cat + b4*JASPfile$rs8139770_cat + b5*JASPfile$rs4821986_cat
#                 + b6*JASPfile$rs4821995_cat + b7*JASPfile$rs114309992_cat + b8*JASPfile$rs20554_cat
#                 + b9*JASPfile$rs1033611_cat + b10*JASPfile$rs2413639_cat + b11*JASPfile$rs139533_cat
#                 + b12*JASPfile$rs11090060_cat + b13*JASPfile$rs1983554_cat + b14*JASPfile$rs9607848_cat
#                 + b15*JASPfile$rs9607850_cat + b16*JASPfile$rs2982055_cat + b17*JASPfile$rs811668_cat
#                 + b18*JASPfile$PC1 + b19*JASPfile$PC2 + b20*JASPfile$PC3 + b21*JASPfile$PC4 + b22*JASPfile$PC5
#                 + b23*JASPfile$mother_age_recruitment + b24*JASPfile$mother_income_cat
#                 + b25*JASPfile$household_income_cat + b26*JASPfile$accommodation_cat
#                 + b27*JASPfile$mother_ethnicity_cat + b28*JASPfile$mother_highest_education_cat
#                 + b29*JASPfile$child.s_sex_binary)

predictions <- (b0+
                  b1*JASPfile$rs133335_cat+
                  b2*JASPfile$rs133344_cat+
                  b3*JASPfile$rs73885718_cat+
                  b4*JASPfile$rs5758165_cat+
                  b5*JASPfile$rs5758166_cat+
                  b6*JASPfile$rs114309992_cat+
                  b7*JASPfile$rs76550409_cat+
                  b8*JASPfile$rs5751117_cat+
                  b9*JASPfile$rs5758487_cat+
                  b10*JASPfile$rs4269017_cat+
                  b11*JASPfile$rs762995_cat+
                  b12*JASPfile$rs134906_cat+
                  b13*JASPfile$rs6002674_cat+
                  b14*JASPfile$rs139533_cat+
                  b15*JASPfile$rs5758689_cat+
                  b16*JASPfile$PC1.x+
                  b17*JASPfile$PC2.x+
                  b18*JASPfile$PC3.x+
                  b19*JASPfile$PC4+
                  b20*JASPfile$PC5+
                  b21*JASPfile$mother_age_recruitment+
                  b22*JASPfile$mother_income_cat+
                  b23*JASPfile$household_income_cat+
                  b24*JASPfile$accommodation_cat+
                  b25*JASPfile$mother_highest_education_cat+
                  b26*JASPfile$mother_chinesevmalay+
                  b27*JASPfile$mother_chinesevindian+
                  b28*JASPfile$child.s_sex_binary)

predictions

# change cgxxxx according to probe number ID
actual <- JASPfile$cg15597984  

# calculate the RMSE
rmse <- sqrt(mean((predictions - actual)^2))
rmse

# calculate the R-squared
ss_res <- sum((predictions - actual)^2)
ss_tot <- sum((actual - mean(actual))^2)
r_squared <- 1 - ss_res/ss_tot
r_squared

#get actual cpg values
cat( paste( JASPfile$cg15597984, collapse='\n' ) )

#get predictions in vertical form for copy & paste
cat( paste( predictions, collapse='\n' ) )
