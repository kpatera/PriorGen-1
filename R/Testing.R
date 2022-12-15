test_mat<-matrix(NA,nrow = 100, ncol = 5)
count<-1
for(i in seq(from =0, to=0.99, length=100)){
  test_mat[count,]<-c(i,findbetaqq(percentile.value1 = 0.4, percentile.value2 = 0.9,
           percentile1 = i, percentile2 = 0.99))
  print(i);count=count+1
}

-1.693147
