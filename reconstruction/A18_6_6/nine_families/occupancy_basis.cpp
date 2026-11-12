// Exact coefficient basis for SPEC01 section 10.
// Does not read model labels or entries. Inputs: independently generated word
// orbit representatives. Counts actual coordinate subsets T without sampling.
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>
using i64 = std::int64_t;
struct Orbit { int v,k; std::array<unsigned,4> words{}; };
using Basis = std::vector<std::array<i64,3>>;
static i64 choose(int n,int k) {
    if(k<0 || k>n) return 0;
    i64 c=1; for(int i=1;i<=k;++i) c=c*(n-i+1)/i; return c;
}
static int power7(int k) {int z=1; while(k--) z*=7; return z;}
static Basis difference_basis(int code,int k,int s,bool complement,int lo,int rows) {
    std::array<int,4> classes{};
    for(int i=0;i<k;++i) {
        int a=code%7; code/=7;
        classes[i]=complement?s-a:a;
        assert(classes[i]>=lo && classes[i]-lo+1<rows);
    }
    Basis out(rows*rows,{0,0,0});
    for(int mask=0;mask<(1<<k);++mask) {
        std::array<int,8> f{}; f[0]=1;
        int n=0;
        for(int z=0;z<k;++z) if(mask&(1<<z)) {
            ++f[classes[z]-lo+1];
            n += classes[z]==s;
        }
        int sign=((k-__builtin_popcount((unsigned)mask))%2)?-1:1;
        std::array<int,3> np{1,n,n*n};
        for(int i=0;i<rows;++i) for(int j=i;j<rows;++j)
            for(int d=0;d<3;++d)
                out[i*rows+j][d]+=(i64)sign*f[i]*f[j]*np[d];
    }
    // Degree of N^d f_i f_j is at most d+2, including Boolean reduction.
    for(const auto& v:out) for(int d=0;d<3;++d)
        if(k>d+2) assert(v[d]==0);
    return out;
}
int main(int argc,char** argv) {
    if(argc!=4) {std::cerr<<"usage: occupancy_basis INPUT OUTPUT TIMINGS\n";return 2;}
    std::ifstream in(argv[1]); std::ofstream out(argv[2]), times(argv[3]);
    int nvars; in>>nvars; assert(nvars==941);
    std::vector<Orbit> orbits;
    for(int z=0;z<nvars;++z) {
        Orbit o; in>>o.v>>o.k; assert(o.v==z+1 && o.k>=1 && o.k<=4);
        for(int i=0;i<o.k;++i) {in>>o.words[i];assert(__builtin_popcount(o.words[i])==6);}
        orbits.push_back(o);
    }
    assert(in.good() || in.eof());
    const std::array<int,5> cores{2,3,7,8,9};
    for(int s:cores) {
        const auto begin=std::chrono::steady_clock::now();
        bool complement=s>=7;
        int weight=complement?12:6;
        int lo=std::max(0,weight-(18-s)), hi=std::min(s,weight),rows=hi-lo+2;
        std::vector<unsigned> ts;
        for(unsigned t=0;t<(1u<<18);++t) if(__builtin_popcount(t)==s) ts.push_back(t);
        assert((i64)ts.size()==choose(18,s));
        std::array<std::map<int,Basis>,5> memo;
        // Empty support: constant component, basis polynomial 1 only.
        out<<s<<' '<<int(complement)<<" 0 0 0 "<<ts.size()<<" 0 0\n";
        for(const Orbit& orb:orbits) {
            std::vector<i64> histogram(power7(orb.k),0);
            for(unsigned t:ts) {
                std::array<int,4> a{};
                for(int w=0;w<orb.k;++w) a[w]=__builtin_popcount(t&orb.words[w]);
                std::sort(a.begin(),a.begin()+orb.k);
                int code=0;for(int w=orb.k-1;w>=0;--w) code=code*7+a[w];
                ++histogram[code];
            }
            assert(std::accumulate(histogram.begin(),histogram.end(),i64(0))==(i64)ts.size());
            Basis totals(rows*rows,{0,0,0});
            for(int code=0;code<(int)histogram.size();++code) if(histogram[code]) {
                auto it=memo[orb.k].find(code);
                if(it==memo[orb.k].end())
                    it=memo[orb.k].emplace(code,difference_basis(code,orb.k,s,complement,lo,rows)).first;
                const auto& part=it->second;
                for(int ij=0;ij<rows*rows;++ij) for(int d=0;d<3;++d)
                    totals[ij][d]+=histogram[code]*part[ij][d];
            }
            for(int i=0;i<rows;++i) for(int j=i;j<rows;++j) {
                const auto& c=totals[i*rows+j];
                if(c[0] || c[1] || c[2])
                    out<<s<<' '<<int(complement)<<' '<<orb.v<<' '<<i<<' '<<j<<' '
                       <<c[0]<<' '<<c[1]<<' '<<c[2]<<'\n';
            }
        }
        auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
        times<<s<<' '<<elapsed<<' '<<ts.size()<<' '<<(i64)ts.size()*nvars<<'\n';
        std::cerr<<"core="<<s<<" complete, coordinate_subsets="<<ts.size()
                 <<", orbit_subset_pairs="<<(i64)ts.size()*nvars<<", seconds="<<elapsed<<'\n';
    }
    return 0;
}
