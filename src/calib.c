#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
typedef uint64_t u64; typedef uint32_t u32; typedef uint8_t u8;
static u32 *primes_upto(u64 limit,size_t*c){size_t n=limit+1;u8*comp=calloc(n,1);size_t cap=limit/10+64;u32*pr=malloc(sizeof(u32)*cap);size_t np=0;
 for(u64 i=2;i<=limit;i++){if(comp[i])continue;if(np<cap)pr[np++]=(u32)i;if(i*i<=limit)for(u64 j=i*i;j<=limit;j+=i)comp[j]=1;}free(comp);*c=np;return pr;}
static double SS(u64 n,u32*base,size_t nb){double s=1;u64 m=n;while(m%2==0)m/=2;for(size_t i=0;i<nb;i++){u64 p=base[i];if(p==2)continue;if(p*p>m)break;if(m%p==0){s*=(double)(p-1)/(p-2);while(m%p==0)m/=p;}}if(m>2)s*=(double)(m-1)/(m-2);return s;}
static double C2=0.66016181584686957392781211001455577843262336;
int main(int argc,char**argv){
 u64 X=argc>=2?strtoull(argv[1],0,10):1000000;
 size_t np=0;u32*pr=primes_upto(X,&np);u8*ip=calloc(X+1,1);for(size_t i=0;i<np;i++)ip[pr[i]]=1;
 size_t nb=0;u32*base=primes_upto(sqrt((double)X)+2,&nb);
 printf("%-12s %-10s %-12s %-12s %-10s\n","n","R(n)","SS(n)","n*S/log^2","R/(n*S/log^2)");
 for(u64 n=X/2;n<=X;n+=X/8){
   u64 R=0; for(size_t i=0;i<np;i++){u64 p=pr[i];if(p>n/2)break;if(ip[n-p])R++;}
   double L=log((double)n); double s=SS(n,base,nb); double den=(double)n*s/(L*L);
   printf("%-12llu %-10llu %-12.4f %-12.1f %-10.4f\n",(unsigned long long)n,(unsigned long long)R,s,den,den>0?(double)R/den:0);
 }
 /* estimate the limiting per-n constant: mean of R(n)*log^2 n/(n*S(n)) over a band */
 double num=0,den2=0;
 for(u64 n=X/2;n<=X;n+=2){u64 R=0;for(size_t i=0;i<np;i++){u64 p=pr[i];if(p>n/2)break;if(ip[n-p])R++;}double L=log((double)n);double s=SS(n,base,nb);num+=(double)R*L*L;den2+=(double)n*s;}
 printf("\nband-averaged per-n constant K = %.5f   (2*C2 = %.5f, C2 = %.5f)\n",num/den2,2*C2,C2);
 printf("=> unordered count matches HL with coefficient K; compare to C2=%.5f\n",C2);
 return 0;}
