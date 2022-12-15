options(scipen = 9999)
set.seed(43241)
First1<-rbinom(n = 1, size = 500, prob = 0.1)

delta=seq(0.01,2,0.01)
std=seq(0.01,2,0.01)
combs=expand.grid(delta,std)
test1=NULL;test1$n=0;test1$delta=0
out1=matrix(NA,ncol = 3,nrow = dim(combs)[1])
for(i in 1:dim(combs)[1]){
    test1=power.t.test(delta = combs[i,1],sd = combs[i,2],sig.level = 0.05,power = 0.8)
    out1[i,]=c(test1$n,test1$delta,test1$sd)
}
out1<-data.frame(out1)
names(out1)<-c("n","m","s")
out_sel=out1[(out1$n>=First1-5 & out1$n<=First1+5 & out1$m>0.2  & out1$m<1.2 & out1$s>0.5 & out1$s<1),]
dim(out_sel)
apply(out_sel,2,mean)

sel_stud<-out_sel[sample(dim(out_sel)[1],1),]

# Option 1 # Larger delta
test_opt1=power.t.test(delta = 0.6,sd = 1,sig.level = 0.05,power = 0.8)
sel_stud$n-test_opt1$n # 4 animals saved.

# Option 2 # Smaller delta
test_opt2=power.t.test(delta = 0.3,sd = 1,sig.level = 0.05,power = 0.8)
sel_stud$n-test_opt2$n # 126 more animals to prove a biologically significant effect
# 48 animals saved.

# Option 3 # Larger delta, smaller st.dev
test_opt3=power.t.test(delta = 0.6,sd = 0.5,sig.level = 0.05,power = 0.8)
sel_stud$n-test_opt3$n # 48 animals saved if bibliography and expectations are properly weighted.


delta=seq(2,0.01,-0.01)
std=seq(2,0.01,-0.01)
combs=expand.grid(delta,std);i=1
test2=NULL;test2$n=0;test2$delta=0
out2=matrix(NA,nrow = 3,ncol = dim(combs)[1])
while(!(test2$n>First1-5 & test2$n<First1+5 & test2$delta>=0.2)){
  test2=power.t.test(delta = combs[i,1],sd = combs[i,2],sig.level = 0.05,power = 0.8)
  out1[,i]=c(test2$n,test2$delta,test2$sd)
  }
out2