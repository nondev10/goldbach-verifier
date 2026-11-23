#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
typedef uint64_t u64; typedef uint32_t u32; typedef uint8_t u8;
static u32 *primes_upto(u64 limit,size_t*c){size_t n=limit+1;u8*comp=calloc(n,1);u32*pr=malloc(sizeof(u32)*(limit/8+64));size_t np=0;
 for(u64 i=2;i<=limit;i++){if(comp[i])continue;pr[np++]=(u32)i;if(i*i<=limit)for(u64 j=i*i;j<=limit;j+=i)comp[j]=1;}free(comp);*c=np;return pr;}
static double SS(u64 n,u32*base,size_t nb){double s=1;u64 m=n;while(m%2==0)m/=2;for(size_t i=0;i<nb;i++){u64 p=base[i];if(p==2)continue;if(p*p>m)break;if(m%p==0){s*=(double)(p-1)/(p-2);while(m%p==0)m/=p;}}if(m>2)s*=(double)(m-1)/(m-2);return s;}
int main(int argc,char**argv){
 u64 x=argc>=2?strtoull(argv[1],0,10):200000;
 size_t np=0;u32*pr=primes_upto(x,&np);u8*ip=calloc(x+1,1);for(size_t i=0;i<np;i++)ip[pr[i]]=1;
 size_t nb=0;u32*base=primes_upto(sqrt((double)x)+2,&nb);
 double sumSS=0;u64 cnt=0;
 for(u64 n=4;n<=x;n+=2){sumSS+=SS(n,base,nb);cnt++;}
 printf("x=%llu  mean SS over even n = %.6f   (C2=0.660162)\n",(unsigned long long)x,sumSS/cnt);
 /* sum of n*SS/log^2 n */
 double A=0; for(u64 n=4;n<=x;n+=2){double L=log((double)n);A+=(double)n*SS(n,base,nb)/(L*L);}
 /* exact total representations */
 double R=0; for(u64 n=4;n<=x;n+=2){u64 c=0;for(size_t i=0;i<np;i++){u64 p=pr[i];if(p>n/2)break;if(ip[n-p])c++;}R+=c;}
 double L=log((double)x);
 printf("sum_R = %.0f\n",R);
 printf("A = sum n*SS/log^2 n = %.4e\n",A);
 printf("implied K = sum_R / A = %.6f   (theory: 2*C2 = %.6f)\n",R/A,2*0.66016181584686957392781211001455577843262336);
 printf("aggregate theory: C2*x^2/(2 log^2 x) = %.4e  ratio=%.4f\n",0.66016181584686957392781211001455577843262336*(double)x*x/(2*L*L), R/(0.66016181584686957392781211001455577843262336*(double)x*x/(2*L*L)));
 return 0;}
