// Independent D0 reconstruction from SPEC01 sections 1--3.
// This program reads only tableaux and an independently enumerated orbit map.
// It never reads model entries, block sizes, or positive_scale_removed.
// All mathematical operations are checked signed 128-bit integer operations.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>
using U = unsigned __int128;
using Z = __int128;
struct Hash {
  size_t operator()(U x) const {
    uint64_t a=uint64_t(x), b=uint64_t(x>>64);
    a ^= b*0x9e3779b97f4a7c15ULL;
    a ^= a>>30; a*=0xbf58476d1ce4e5b9ULL;
    a ^= a>>27; a*=0x94d049bb133111ebULL; return a^(a>>31);
  }
};
using Poly = std::unordered_map<U,Z,Hash>;
Z add(Z a,Z b) { Z r; if(__builtin_add_overflow(a,b,&r)) throw std::overflow_error("integer addition"); return r; }
Z mul(Z a,Z b) { Z r; if(__builtin_mul_overflow(a,b,&r)) throw std::overflow_error("integer multiplication"); return r; }
std::string str(Z a) { if(!a)return "0"; bool neg=a<0; if(neg)a=-a; std::string s; while(a){s.push_back(char('0'+a%10));a/=10;} if(neg)s.push_back('-');std::reverse(s.begin(),s.end());return s; }
void inc(Poly &p,U h,Z c) { if(!c)return; auto it=p.find(h); if(it==p.end())p.emplace(h,c); else it->second=add(it->second,c); }
void clean(Poly &p) { for(auto it=p.begin();it!=p.end();)if(!it->second)it=p.erase(it);else ++it; }
int expn(U h,int a){return int((h>>(5*a))&31);}
U unit(int a){return U(1)<<(5*a);}
long long choose[19][19];
size_t max_terms=0, final_terms=0, allowed_terms=0, pairs=0, lowering_calls=0;
Poly product(const Poly&a,const Poly&b){ Poly out; for(const auto&[x,c]:a)for(const auto&[y,d]:b)inc(out,x+y,mul(c,d));clean(out);max_terms=std::max(max_terms,out.size());return out; }
Poly highest(const std::array<int,4>&lambda){
  Poly out{{0,1}};
  for(int l=1;l<=4;l++){
    int power=lambda[l-1]-(l==4?0:lambda[l]);
    if(!power)continue;
    Poly det;std::vector<int> p(l);std::iota(p.begin(),p.end(),0);
    int factorial=1;for(int i=2;i<=l;i++)factorial*=i;
    do{int inv=0;U h=0;for(int a=0;a<l;a++){h+=unit(4*a+p[a]);for(int b=a+1;b<l;b++)inv+=(p[a]>p[b]);}det[h]=(inv%2?-factorial:factorial);}while(std::next_permutation(p.begin(),p.end()));
    while(power--)out=product(out,det);
  }
  return out;
}
// Divided power E_ij^e/e!: distribute e moves among four source cells.
// The coefficient is a product of ordinary binomial coefficients.
Poly lower(const Poly&in,int source,int target,int e,bool right){
  if(!e)return in;
  ++lowering_calls;
  Poly out;out.reserve(in.size()*2+1);
  int src[4],dst[4];U su[4],du[4];
  for(int b=0;b<4;b++){src[b]=right?4*b+source:4*source+b;dst[b]=right?4*b+target:4*target+b;su[b]=unit(src[b]);du[b]=unit(dst[b]);}
  for(const auto&[h,c]:in){
    int a[4];for(int b=0;b<4;b++)a[b]=expn(h,src[b]);
    for(int q0=std::max(0,e-a[1]-a[2]-a[3]);q0<=std::min(a[0],e);q0++)
    for(int q1=std::max(0,e-q0-a[2]-a[3]);q1<=std::min(a[1],e-q0);q1++)
    for(int q2=std::max(0,e-q0-q1-a[3]);q2<=std::min(a[2],e-q0-q1);q2++){
      int q3=e-q0-q1-q2;if(q3<0||q3>a[3])continue;
      int q[4]={q0,q1,q2,q3};U key=h;Z mult=1;
      for(int b=0;b<4;b++){key-=su[b]*q[b];key+=du[b]*q[b];mult=mul(mult,choose[a[b]][q[b]]);}
      inc(out,key,mul(c,mult));
    }
  }
  clean(out);max_terms=std::max(max_terms,out.size());return out;
}
using Op=std::tuple<int,int,int>;
struct Node {std::map<Op,int> children;std::vector<int> leaves;int maxleaf=-1;};
int main(int argc,char**argv){
 try{
  if(argc!=5)throw std::runtime_error("usage: d0_engine D0_LABELS_TXT ORBIT_LOOKUP_TSV BLOCK OUT_PREFIX");
  auto started=std::chrono::steady_clock::now();auto cpu=std::clock();
  for(int n=0;n<=18;n++){choose[n][0]=choose[n][n]=1;for(int k=1;k<n;k++)choose[n][k]=choose[n-1][k-1]+choose[n-1][k];}
  int wanted=std::stoi(argv[3]);std::ifstream labels(argv[1]);if(!labels)throw std::runtime_error("labels file");
  std::array<int,4> lambda;std::vector<std::vector<Op>> ops;std::vector<std::array<int,4>> contents;
  int bid,n;
  while(labels>>bid>>n){
    std::array<int,4> shape;for(int &a:shape)labels>>a;
    std::vector<std::vector<Op>> blockops;std::vector<std::array<int,4>> blockcontents;
    for(int t=0;t<n;t++){
      int counts[4][4]{};std::array<int,4> content{};
      for(int row=0;row<4;row++)for(int i=0;i<shape[row];i++){int a;labels>>a;counts[row][a]++;content[a]++;}
      std::vector<Op> v;for(int j=3;j>=0;j--)for(int i=3;i>j;i--)if(counts[j][i])v.emplace_back(j,i,counts[j][i]);
      blockops.push_back(v);blockcontents.push_back(content);
    }
    if(bid==wanted){lambda=shape;ops=std::move(blockops);contents=std::move(blockcontents);break;}
  }
  if(ops.empty())throw std::runtime_error("block not found");
  std::unordered_map<U,int,Hash> lookup;std::ifstream lut(argv[2]);
  if(!lut)throw std::runtime_error("orbit lookup file");
  while(true){int a,v;U key=0;if(!(lut>>a))break;key=U(a);for(int i=1;i<16;i++){lut>>a;key+=U(a)<<(5*i);}lut>>v;if(v<1||v>941)throw std::runtime_error("lookup variable");lookup.emplace(key,v);}
  std::string prefix=argv[4];std::ofstream out(prefix+".raw.tsv");if(!out)throw std::runtime_error("output file");
  int offset=(lambda[0]==18);n=ops.size();
  if(offset){
    out<<"0 0 0 18564\n";
    Z fact=1;for(int k=1;k<=18;k++)fact*=k;
    for(int i=0;i<n;i++){
      Z denom=1;for(int a:contents[i])for(int k=1;k<=a;k++)denom*=k;
      if(fact%denom)throw std::runtime_error("multinomial not integral");
      U hist=0;for(int a=0;a<4;a++)hist+=U(contents[i][a])<<(5*(4*a+a));
      auto f=lookup.find(hist);if(f==lookup.end())throw std::runtime_error("single-row content lookup");
      out<<f->second<<" 0 "<<i+1<<" "<<str(fact/denom)<<"\n";
    }
  }
  Poly h=highest(lambda);
  std::vector<Node> trie(1);
  for(int j=0;j<n;j++){
    int at=0;trie[at].maxleaf=std::max(trie[at].maxleaf,j);
    for(auto op:ops[j]){
      auto f=trie[at].children.find(op);int next;
      if(f==trie[at].children.end()){next=trie.size();trie[at].children[op]=next;trie.emplace_back();}else next=f->second;
      at=next;trie[at].maxleaf=std::max(trie[at].maxleaf,j);
    }
    trie[at].leaves.push_back(j);
  }
  for(int i=0;i<n;i++){
    Poly left=h;for(auto [s,t,e]:ops[i])left=lower(left,s,t,e,false);
    std::function<void(int,const Poly&)> visit=[&](int at,const Poly&p){
      for(int j:trie[at].leaves)if(j>=i){
        std::array<Z,942> coeff{};pairs++;final_terms+=p.size();
        for(const auto&[key,c]:p){auto f=lookup.find(key);if(f!=lookup.end()){coeff[f->second]=add(coeff[f->second],c);allowed_terms++;}}
        for(int v=1;v<=941;v++)if(coeff[v])out<<v<<" "<<i+offset<<" "<<j+offset<<" "<<str(coeff[v])<<"\n";
      }
      for(const auto&[op,next]:trie[at].children)if(trie[next].maxleaf>=i){auto [s,t,e]=op;Poly q=lower(p,s,t,e,true);visit(next,q);}
    };
    visit(0,left);
    if(i==0||(i+1)%10==0||i+1==n){double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();std::cout<<"D0 "<<wanted<<" left="<<i+1<<"/"<<n<<" elapsed="<<sec<<" max_terms="<<max_terms<<std::endl;}
  }
  out.close();double wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();double cpu_sec=double(std::clock()-cpu)/CLOCKS_PER_SEC;
  std::ofstream metrics(prefix+".metrics.json");
  metrics<<"{\"block\":"<<wanted<<",\"tableau_columns\":"<<n<<",\"full_columns\":"<<n+offset<<",\"tableau_pairs\":"<<pairs<<",\"max_polynomial_terms\":"<<max_terms<<",\"total_final_terms\":"<<final_terms<<",\"allowed_final_terms\":"<<allowed_terms<<",\"lowering_calls\":"<<lowering_calls<<",\"wall_seconds\":"<<wall<<",\"cpu_seconds\":"<<cpu_sec<<",\"integer_overflow\":false}\n";
  return 0;
 }catch(const std::exception&e){std::cerr<<"ERROR "<<e.what()<<std::endl;return 1;}
}
