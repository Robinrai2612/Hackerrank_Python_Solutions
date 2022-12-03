#include<bits/stdc++.h>
using namespace std;

#define REP(i,a,b) for(i=a;i<b;i++)
#define rep(i,n) REP(i,0,n)

#define mygc(c) (c)=getchar_unlocked()
#define mypc(c) putchar_unlocked(c)

#define ll long long
#define ull unsigned ll

void reader(int *x){int k,m=0;*x=0;for(;;){mygc(k);if(k=='-'){m=1;break;}if('0'<=k&&k<='9'){*x=k-'0';break;}}for(;;){mygc(k);if(k<'0'||k>'9')break;*x=(*x)*10+k-'0';}if(m)(*x)=-(*x);}
void reader(ll *x){int k,m=0;*x=0;for(;;){mygc(k);if(k=='-'){m=1;break;}if('0'<=k&&k<='9'){*x=k-'0';break;}}for(;;){mygc(k);if(k<'0'||k>'9')break;*x=(*x)*10+k-'0';}if(m)(*x)=-(*x);}
void reader(double *x){scanf("%lf",x);}
int reader(char c[]){int i,s=0;for(;;){mygc(i);if(i!=' '&&i!='\n'&&i!='\r'&&i!='\t'&&i!=EOF) break;}c[s++]=i;for(;;){mygc(i);if(i==' '||i=='\n'||i=='\r'||i=='\t'||i==EOF) break;c[s++]=i;}c[s]='\0';return s;}
template <class T, class S> void reader(T *x, S *y){reader(x);reader(y);}
template <class T, class S, class U> void reader(T *x, S *y, U *z){reader(x);reader(y);reader(z);}
template <class T, class S, class U, class V> void reader(T *x, S *y, U *z, V *w){reader(x);reader(y);reader(z);reader(w);}

void writer(int x, char c){int s=0,m=0;char f[10];if(x<0)m=1,x=-x;while(x)f[s++]=x%10,x/=10;if(!s)f[s++]=0;if(m)mypc('-');while(s--)mypc(f[s]+'0');mypc(c);}
void writer(ll x, char c){int s=0,m=0;char f[20];if(x<0)m=1,x=-x;while(x)f[s++]=x%10,x/=10;if(!s)f[s++]=0;if(m)mypc('-');while(s--)mypc(f[s]+'0');mypc(c);}
void writer(double x, char c){printf("%.15f",x);mypc(c);}
void writer(const char c[]){int i;for(i=0;c[i]!='\0';i++)mypc(c[i]);}
void writer(const char x[], char c){int i;for(i=0;x[i]!='\0';i++)mypc(x[i]);mypc(c);}
template<class T> void writerLn(T x){writer(x,'\n');}
template<class T, class S> void writerLn(T x, S y){writer(x,' ');writer(y,'\n');}
template<class T, class S, class U> void writerLn(T x, S y, U z){writer(x,' ');writer(y,' ');writer(z,'\n');}
template<class T> void writerArr(T x[], int n){int i;if(!n){mypc('\n');return;}rep(i,n-1)writer(x[i],' ');writer(x[n-1],'\n');}

template<class T> void sort(int N, T a[], void *mem = NULL){sort(a,a+N);}
template<class T1, class T2> void sort(int N, T1 a[], T2 b[], void *mem){int i;pair<T1,T2> *r=(pair<T1, T2>*)mem;rep(i,N)r[i].first=a[i],r[i].second=b[i];sort(r,r+N);rep(i,N)a[i]=r[i].first,b[i]=r[i].second;}
template<class T1, class T2, class T3> void sort(int N, T1 a[], T2 b[], T3 c[], void *mem){int i;pair<T1,pair<T2,T3> > *r=(pair<T1,pair<T2,T3> >*)mem;rep(i,N)r[i].first=a[i],r[i].second.first=b[i],r[i].second.second=c[i];sort(r,r+N);rep(i,N)a[i]=r[i].first,b[i]=r[i].second.first,c[i]=r[i].second.second;}
template<class T1, class T2, class T3, class T4> void sort(int N, T1 a[], T2 b[], T3 c[], T4 d[], void *mem){int i;pair<pair<T1,T2>,pair<T3,T4> > *r=(pair<pair<T1,T2>,pair<T3,T4> >*)mem;rep(i,N)r[i].first.first=a[i],r[i].first.second=b[i],r[i].second.first=c[i],r[i].second.second=d[i];sort(r,r+N);rep(i,N)a[i]=r[i].first.first,b[i]=r[i].first.second,c[i]=r[i].second.first,d[i]=r[i].second.second;}


char memarr[77000000]; void *mem = memarr;
#define MD 1000000007

template<class T>
struct rollingHash64{
  int len;
  T *data;
  ull *sum, *rev, *pw;
  ull mul;

  ull getinv(ull a){
    ull t,s=a,u=0,v=1,e;
    e = numeric_limits<ull>::max() / s;
    t -= e * s;
    u -= e * v;
    swap(t,s);
    swap(u,v);
    while(s){
      e=t/s;
      t-=e*s;
      u-=e*v;
      swap(t,s);
      swap(u,v);
    }
    return u;
  }

  void* init(int n, T *arr, ull m = 0, void *mem = NULL){
    int i; ull v;

    mul = m;
    if(mul==0) mul = 2*(rand()%1000000000) + 1000000001ULL;

    len = n;
    data = arr;
    if(mem == NULL){
      pw = (ull*)malloc(sizeof(ull)*(2*len+1));
      sum = (ull*)malloc(sizeof(ull)*(len+1));
      rev = (ull*)malloc(sizeof(ull)*(len+1));
    } else {
      pw = (ull*)mem;
      sum = pw + 2*len + 1;
      rev = sum + len + 1;
      mem = rev + len + 1;
    }

    v = getinv(mul);
    pw = pw + len;
    pw[0] = 1;
    rep(i,len) pw[ i+1] = pw[ i] * mul;
    rep(i,len) pw[-i-1] = pw[-i] * v;

    sum[0] = 0;
    rep(i,len) sum[i+1] = sum[i] + (ull)data[i] * pw[i];

    rev[len] = 0;
    for(i=len-1;i>=0;i--) rev[i] = rev[i+1] + (ull)data[i] * pw[len-i-1];

    return mem;
  }

  ull get(int a, int b, int off=0){
    ull res;
    
    if(a <= b){
      res = (sum[b+1] - sum[a]) * pw[-a+off] + (b-a+1);
    } else {
      res = (rev[b] - rev[a+1]) * pw[-(len-1-a)+off] + (a-b+1);
    }

    return res;
  }
};

template<class T>
void manacher(int n, T arr[], int res[]) {
  int i, j, k;
  for(i=0,j=0; i<2*n; i+=k, j=max(j-k,0)) {
    while(i-j >= 0 && i+j+1 < 2*n && arr[(i-j)/2] == arr[(i+j+1)/2]) ++j;
    res[i] = j;
    for(k=1; i-k >= 0 && res[i]-k >= 0 && res[i-k] != res[i]-k; ++k)
      res[i+k] = min(res[i-k], res[i]-k);
  }
}


template<class T>
struct lazySegtreeMinVal{
  int N, logN;
  T *data;

  T *fixval; char *fixed;
  T *addval;

  void malloc(int maxN){
    int i;
    for(i=1;i<maxN;i*=2);
    
    data = (T*)std::malloc(sizeof(T)*2*i);
    fixval = (T*)std::malloc(sizeof(T)*i);
    addval = (T*)std::malloc(sizeof(T)*i);
    fixed = (char*)std::malloc(sizeof(char)*i);
  }

  T& operator[](int i){
    return data[N+i];
  }

  void setN(int n, int zerofill = 1){
    int i;
    for(i=1,logN=0;i<n;i*=2,logN++);
    N = i;
    if(zerofill) rep(i,N) data[N+i] = 0;
  }

  void build(void){
    int i;
    for(i=N-1;i;i--) data[i] = min(data[2*i],data[2*i+1]);
    REP(i,1,N) fixed[i] = 0;
    REP(i,1,N) addval[i] = 0;
  }

  inline void push_one(int a, int sz){
    if(fixed[a]){
      if(sz > 1){
        fixed[a*2] = fixed[a*2+1] = 1;
        fixval[a*2] = fixval[a*2+1] = fixval[a];
        data[a*2] = data[a*2+1] = fixval[a];
      } else {
        data[a*2] = data[a*2+1] = fixval[a];
      }
      fixed[a] = 0;
      addval[a] = 0;
      return;
    }
    if(addval[a] != 0){
      if(sz > 1){
        if(fixed[a*2]) fixval[a*2] += addval[a];
        else           addval[a*2] += addval[a];
        if(fixed[a*2+1]) fixval[a*2+1] += addval[a];
        else             addval[a*2+1] += addval[a];
        data[a*2] += addval[a];
        data[a*2+1] += addval[a];
      } else {
        data[a*2] += addval[a];
        data[a*2+1] += addval[a];
      }
      addval[a] = 0;
      return;
    }
  }

  inline void push(int a){
    int i, aa;
    for(i=logN;i;i--){
      aa = a>>i;
      push_one(aa, 1<<(i-1));
    }
  }

  inline void build(int a){
    while(a > 1){
      a /= 2;
      if(fixed[a]){
        data[a] = fixval[a];
      } else {
        data[a] = min(data[a*2], data[a*2+1]);
        if(addval[a] != 0) data[a] += addval[a];
      }
    }
  }

  inline void change(int a, int b, T val){
    int aa, bb;
    if(a >= b) return;

    aa = (a += N);
    bb = (b += N);
    push(a); push(b-1);

    if(a%2) data[a++] = val;
    if(b%2) data[--b] = val;
    a /= 2;
    b /= 2;

    while(a < b){
      if(a%2) fixed[a]=1, fixval[a]=val, data[a++] = val;
      if(b%2) fixed[--b]=1, fixval[b]=val, data[b] = val;
      a /= 2;
      b /= 2;
    }

    build(aa);
    build(bb-1);
  }

  inline void add(int a, int b, T val){
    int sz = 1, aa, bb;
    if(a >= b) return;

    aa = (a += N);
    bb = (b += N);
    push(a); push(b-1);

    if(a%2) data[a++] += val;
    if(b%2) data[--b] += val;
    a /= 2;
    b /= 2;

    while(a < b){
      sz *= 2;
      if(a%2){
        if(fixed[a]) fixval[a] += val; else addval[a] += val;
        data[a++] += val;
      }
      if(b%2){
        b--;
        if(fixed[b]) fixval[b] += val; else addval[b] += val;
        data[b] += val;
      }
      a /= 2;
      b /= 2;
    }

    build(aa);
    build(bb-1);
  }

  inline T getMinVal(int a, int b){
    T res;
    int sz = 1;
    
    a += N;
    b += N;
    push(a); push(b-1);

    res = std::numeric_limits<T>::max();
    while(a < b){
      if(a%2) res = min(res, data[a++]);
      if(b%2) res = min(res, data[--b]);
      a /= 2;
      b /= 2;
    }
    return res;
  }
};


int N;
char S[2000000];
int rad[3000000];
int res[1000000];

int ss[3000000], ee[3000000], vv[3000000], nx[1000000];

int get_nx(int i){
  if(nx[i]==-1) return i;
  if(i==N-1) return nx[i] = N;
  return nx[i] = get_nx(nx[i]);
}

int main(){
  int i, j, k, st, ed, m, d;
//  rollingHash64<char> h;
//  lazySegtreeMinVal<int> t;

  reader(&N,S);
//  h.init(N,S);
//  t.malloc(N);
//  t.setN(N);
//  t.build();

  rep(i,N) S[N+i] = S[i];

  manacher(2*N, S, rad);
  rep(i,4*N){
    k = min(N,rad[i]);
    if(i%2==0 && k%2==0) k--;
    if(i%2==1 && k%2==1) k--;
    if(rad[i]==0) continue;
    st = i/2 - (k-1)/2;
    ed = i/2 + k/2;

    m = ed-st+1;

    ss[i] = (st-(N-m)+N+N)%N;
    ee[i] = (st+N+N)%N;
    vv[i] = m;
//    rep(j,N-m+1) res[(st+N-j)%N] = max(res[(st+N-j)%N], m);
  }

  sort(4*N, vv, ss, ee, mem);
  rep(i,N+1) nx[i] = -1;
  for(i=4*N-1;i>=0;i--){
    if(ss[i] <= ee[i]){
      k = ss[i];
      while(k <= ee[i]){
//        writerLn(ss[i],k,ee[i]);
        res[k] = max(res[k], vv[i]);
        if(nx[k]==-1) nx[k] = k+1;
        k = get_nx(k);
      }
    } else {
      k = ss[i];
      while(k < N){
//        writerLn(ss[i],k,N);
        res[k] = max(res[k], vv[i]);
        if(nx[k]==-1) nx[k] = k+1;
        k = get_nx(k);
      }

      k = 0;
      while(k <= ee[i]){
//        writerLn(0,k,ee[i]);
        res[k] = max(res[k], vv[i]);
        if(nx[k]==-1) nx[k] = k+1;
        k = get_nx(k);
      }
    }
  }

  REP(i,1,2*N) res[i%N] = max(res[i%N], res[(i-1)%N]-2);
  for(i=2*N-2;i>=0;i--) res[i%N] = max(res[i%N], res[(i+1)%N]-2);

  rep(i,N) writerLn(res[i]);

  return 0;
}
