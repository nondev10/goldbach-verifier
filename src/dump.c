/* dump p_min(n) for even n <= LIMIT, using an INDEPENDENT simple method
   (plain full sieve + direct scan), to cross-check the fast verifier. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
typedef uint64_t u64; typedef uint32_t u32; typedef uint8_t u8;
int main(int argc,char**argv){
 u64 L=argc>=2?strtoull(argv[1],0,10):200000;
 u8*comp=calloc(L+1,1); u8*isp=calloc(L+1,1);
 u32*pr=malloc(sizeof(u32)*(L/8+64)); size_t np=0;
 for(u64 i=2;i<=L;i++){ if(comp[i]) continue; isp[i]=1; pr[np++]=(u32)i;
   if(i*i<=L) for(u64 j=i*i;j<=L;j+=i) comp[j]=1; }
 for(u64 n=4;n<=L;n+=2){ u64 pm=0;
   for(size_t k=0;k<np;k++){u64 p=pr[k]; if(p>n/2)break; if(isp[n-p]){pm=p;break;}}
   printf("%llu %llu\n",(unsigned long long)n,(unsigned long long)pm);}
 return 0;}
