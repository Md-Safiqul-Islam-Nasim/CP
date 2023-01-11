                           ----My Template ---
/**Bismillahir Rahmanir Rahim.**/
/**
 *    author:  nasimnoob
**/
#include<bits/stdc++.h>
using namespace std;
//----------------------------------------------------------//
typedef long long int                        ll;
typedef unsigned long long int              ull;
typedef long double                          ld;
//-----------------------------------------------------------//
#define      test                       int t; cin>>t;  while(t--)
#define      dot(x)                     fixed << setprecision(x)
#define      Case                       cout << "Case " << ++tt << ": " ;
#define      PI                         acos(-1) // 3.1415926535897932
#define      EPS                        (1e-7)
#define      Nasim_Noob                 ios_base::sync_with_stdio(0);cin.tie(0);cout.tie(0);
#define      Mod                        1000000007

template <typename T> T Gcd(T a, T b) {if (a < 0)return Gcd(-a, b); if (b < 0)return Gcd(a, -b); return (b == 0) ? a : Gcd(b, a % b);} // better than __gcd
template <typename T> T Lcm(T a, T b) {if (a < 0)return Lcm(-a, b); if (b < 0)return Lcm(a, -b); return a * (b / Gcd(a, b));}
ll BigMod(ll p, ll e, ll M) {ll ret = 1; for (; e > 0; e >>= 1) { if (e & 1) ret = (ret * p) % M; p = (p * p) % M;} return (ll)ret;}
ll modInverse(ll a, ll M) {return BigMod(a, M - 2, M);}

// priority_queue<int> Q (greatar to smaller)
// priority_queue<int, vector<int>, greater<int> > Q;   (smaller to greater)

int chessx[] = { -1, -1, 1, 1, -2, -2, 2, 2}; //knight
int chessy[] = { 2, -2, 2, -2, 1, -1, -1, 1}; //knight
int xx[] = { 0, 0, 1, -1};
int yy[] = { -1, 1, 0, 0};
//scanf("%[^\n]%*c", str);

//--------------------------------------Solution------------------------------//
int tt = 0;

void solve()
{
  test{
    
  }
}

int main()
{
#ifdef NASIM_NOOB
  clock_t tStart = clock();
  freopen("input.txt", "r", stdin);
  freopen("output.txt", "w", stdout);
#endif

  Nasim_Noob ;
  solve();


#ifdef NASIM_NOOB
  fprintf(stderr, "\n>> Runtime: %.10fs\n", (double) (clock() - tStart) / CLOCKS_PER_SEC);
#endif

  return 0 ;
}

--------------------------------------------------------------------------------------------
                         ------ Algorithm Part ------
--------------------------------------------------------------------------------------------

                         ------ Number Theory -------

        ---- Prime Number ----
------------------------------------------------------------------

bool mark[1000000];
vector < ll > prime ;
void Sieve(ll num)
{
  memset(mark, true, sizeof(mark));

  for (ll i = 2; i * i <= num; i++)
  {
    if (mark[i] == true)
    {
      for (ll j = i * 2 ; j <= num; j += i)
      {
        mark[j] = false;
      }
    }
  }

  for (int i = 2 ; i <= num ; i++)
  {
    if (mark[i])
      prime.push_back(i);
  }
}
------------------------------------------------------------------
     ----- Bitwise Prime -----
------------------------------------------------------------------
vector < ll > Prime;
const int mx = 2 * 1e8 ;
int prime[mx / 32];
void mark_true(int pos){
    prime[pos / 32] |= (1 << (pos % 32));
}
bool is_prime(int pos){
    return !(prime[pos / 32] & (1 << (pos % 32)));
}
void Sieve(){
    int l = sqrt(mx) + 1;
    Prime.push_back(2);
    mark_true(1);
    mark_true(0);
    for (int i = 4 ; i < mx ; i += 2)mark_true(i);
    for (int i = 3 ; i < mx ; i += 2){
        if(is_prime(i)){
            Prime.push_back(i);
            if(i <= l){
                for (int j = i * i ; j < mx ; j += i * 2){
                    mark_true(j);
                }
            }
        }
    }
}
------------------------------------------------------------------


    ----Range of Prime Number----
------------------------------------------------------------------
bool mark[1000000];
vector < ll > prime1 , prime ;

void Sieve(ll num)
{
  memset(mark, true, sizeof(mark));

  for (ll i = 2; i * i <= num; i++)
  {
    if (mark[i] == true)
    {
      for (ll j = i * 2 ; j <= num; j += i)
      {
        mark[j] = false;
      }
    }
  }

  for (int i = 2 ; i <= num ; i++)
  {
    if (mark[i])
      prime1.push_back(i);
  }
}

void primesInRange(ll low, ll high)
{
  ll limit = floor(sqrt(high)) + 1;
  //vector<int> prime;
  Sieve(limit);

  ll n = high - low + 1;

  bool mark[n + 1];
  memset(mark, false, sizeof(mark));

  for (ll i = 0; i < prime1.size(); i++)
  {
    ll loLim = floor(low / prime1[i]) * prime1[i];
    if (loLim < low)
      loLim += prime1[i];
    if (loLim == prime1[i])
      loLim += prime1[i];

    for (ll j = loLim; j <= high; j += prime1[i])
      if (j != prime1[i])
        mark[j - low] = true;
  }

  for (ll i = low; i <= high; i++)
    if (!mark[i - low])
    {
      if (i != 1)
        prime.push_back(i);
    }
}

-------------------------------------------------------------------
    -----**If I have to face a large number of factorial 
           then i have to use this formula
           firstly i have to find primefactor of that number ------
-------------------------------------------------------------------
vector < ll > primefactor ;

void primeFactors(int n)
{
  while (n % 2 == 0)
  {
    primefactor.push_back(2);
    n = n / 2;
  }

  for (int i = 3; i <= sqrt(n); i = i + 2)
  {
    while (n % i == 0)
    {
      primefactor.push_back(i);
      n = n / i;
    }
  }
  if (n > 2)
    primefactor.push_back(n);
}

for (int i = 0 ; i < primefactor.size() ; i ++ )
{

  ll prime = primefactor[i].first;
  ll factor = primefactor[i].second;
  ll copy = prime ;
  int cnt = 0 ;
  while (n / copy)
  {
    cnt += (n / copy);
    copy *= prime ;
  }

  if (cnt < factor)
  {
    found = false ;
    break;
  }
}
-------------------------------------------------------------------
*** Find first  3  digit of n^k ;
-------------------------------------------------------------------

double p = (double) k * log10(n * 1.0);
p -= (int)p;
double ok = pow(10.0, p);

-------------------------------------------------------------------
          Big Mod
-------------------------------------------------------------------
ll bigmod(ll a , ll p , ll m)
{
  if (p == 0)
    return 1; 
  if (p % 2) 
  {
    return ((a % m) * (bigmod(a , p - 1 , m))) % m;
  }
  else 
  {
    ll c = bigmod (a , p / 2 , m);
    return ((c % m) * (c % m)) % m;
  }
}
-------------------------------------------------------------------

   **Factorial (For Large Number)**

---------------------------------------------------------------------------- -
vector < int > Factorial ;

void multiply(int x)
{
  int carry = 0;
  int len = Factorial.size();

  for (int i = 0; i < len ; i++)
  {
    int prod = Factorial[i] * x + carry;

    Factorial[i] = prod % 10;

    carry  = prod / 10;
  }

  while (carry)
  {
    Factorial.push_back(carry % 10);
    carry = carry / 10;
  }


  for (int i = Factorial.size() - 1 ; i >= 0  ; i --)
    cout << Factorial[i] ;
  cout << endl;

}
---------------------------------------------------------------------------- -
          **Fibonacci (For Large Number)**
---------------------------------------------------------------------------- -
string Fibo[5010] ;

string add(string a , string b)
{
  int sum , carry = 0 ;
  string res ;

  int len  = a.size() ;
  int len1 = b.size();

  if (len > len1)
  {
    reverse(b.begin(), b.end());
    for (int i = 0 ; i < (len - len1) ; i++)
      b += '0' ;

    reverse(b.begin(), b.end());
  }

  else if (len < len1)
  {
    reverse(a.begin(), a.end());
    for (int i = 0 ; i < (len1 - len) ; i++)
      a += '0' ;

    reverse(a.begin(), a.end());
  }

  for (int i = a.size() - 1 ; i >= 0 ; i--)
  {
    sum = a[i] - '0' + b[i] - '0' + carry ;

    res += (sum % 10) + '0' ;
    carry = sum / 10 ;
  }

  while (carry)
  {
    res += carry + '0' ;
    carry /= 10 ;
  }

  reverse(res.begin(), res.end()) ;

  return res;

}
-------------------------------------------------------------------

         **Permutation Without Function**
-------------------------------------------------------------------

void Permuation_WithoutFun(string str , ll n)
{
  int len = str.size() ;

  int cnt = 0;

  vector < char > vec , permut ;

  for (int i = 0 ; i < len ; i++)
  {
    vec.push_back(str[i]);
  }

  /*for (int i = 0 ; i < vec.size() ; i ++)
      cout << vec[i] << " ";
  cout << endl;*/

  ll fact[len] ;

  fact[0] = 1 ;

  sort(vec.begin() , vec.end()) ;

  for (int i = 1 ; i < len ; i++)
  {
    fact[i] = fact[i - 1] * i ;
  }

  /*for (int i = 0 ; i < len ; i++)
      cout << fact[i] << " ";
  cout << endl;*/

  while ( !vec.empty())
  {
    ll x = fact[vec.size() - 1];
    ll y = (n - 1) / x ;

    permut.push_back(vec[y]);
    vec.erase(vec.begin() + y);

    n = n - y * x ;
  }

  for (int i = 0 ; i < permut.size() ; i ++)
  {
    cout << permut[i] ;
  }

  cout << endl;
}

-------------------------------------------------------------------
          NcR
------------------------------------------------------------------
ll printNcR(ll n, ll r)
{
  ll p = 1, k = 1;

  if (n - r < r)
    r = n - r;

  if (r != 0) {
    while (r) {
      p *= n;
      k *= r;

      ll m = __gcd(p, k);
      p /= m;
      k /= m;

      n--;
      r--;
    }
  }

  else
    p = 1;

  return p ;
}

-------------------------------------------------------------------

           //count co-prime of n

-------------------------------------------------------------------
ll phi(ll n)
{
  ll result = n;
  for (ll i = 2; i * i <= n; i++)
  {
    if (n % i == 0)
    {
      while (n % i == 0)
        n /= i;
      result -= result / i;
    }
  }
  if (n > 1)
    result -= result / n;
  return result;
}
-------------------------------------------------------------------
        //count co-prime betwenn 1 to n
-------------------------------------------------------------------
void phi_1_to_n(int n) {
  vector<int> phi(n + 1);
  phi[0] = 0;
  phi[1] = 1;
  for (int i = 2; i <= n; i++)
    phi[i] = i;
  for (int i = 2; i <= n; i++) {
    if (phi[i] == i) {
      for (int j = i; j <= n; j += i)
        phi[j] -= phi[j] / i;
    }
  }
}
-------------------------------------------------------------------
    LCS
-------------------------------------------------------------------
void lcs(string s, string st, int n, int m)
{
  int ar [n + 1][m + 1];

  for (int i = 0; i <= n; i++)
  {
    for (int j = 0; j <= m; j++)
    {
      if (i == 0 || j == 0)
      {
        ar[i][j] = 0;
      }
      else if (s[i - 1] == st[j - 1])
      {
        ar[i][j] = ar[i - 1][j - 1] + 1;
      }
      else
      {
        ar[i][j] = max(ar[i - 1][j], ar[i][j - 1]);
      }
    }
  }

  string str ;

  int i = n, j = m;

  while (i > 0 && j > 0)
  {
    if (s[i - 1] == st[j - 1])
    {
      str += s[i - 1];
      i--;
      j--;
    }
    else if (ar[i - 1][j] > ar[i][j - 1])
      i--;
    else
      j--;
  }

  reverse(str.begin(), str.end());

  cout << str << endl;

  cout << str.size() << endl;
}
-------------------------------------------------------------
         InclusionExclusion
-------------------------------------------------------------
ll InclusionExclusion(ll ar[], ll m , ll n)
{
    ll odd = 0, even = 0 , p = 1;
    ll len = (1 << n);
    for (ll i = 1; i < len ; i ++) {
        p = 1;
        for (ll j = 0; j < n; j++) {
            if (i & (1 << j)) {
                p = Lcm(p , ar[j]);
            }
        }
        if (__builtin_popcount(i) & 1)
            odd += (m / p);
        else
            even += (m / p);
    }
    return odd - even;
}

-------------------------------------------------------------
         Merge Sort
-------------------------------------------------------------

int ar[100005];
ll cnt = 0;
vector < int > lft , rght;
void merge(int l , int r , int mid){
    lft.clear();
    rght.clear();
    int i;
    for(i = l ; i < mid ; i ++){
        lft.push_back(ar[i]);
    }
    for(i = mid ; i < r ; i ++){
        rght.push_back(ar[i]);
    }
    int j = 0 , k ;
    i = 0;
    k = l;
    while(i < lft.size() && j < rght.size()){
          if(lft[i] <= rght[j]){
             ar[k++] = lft[i++];
          }
          else{
              ar[k++] = rght[j++];
          }
    }
    while(i <lft.size()){
         ar[k++] = lft[i++];
    }
    while(j < rght.size()){
         ar[k++] = rght[j++];
    }
}
 
void mergesort(int l, int r){
    if(l + 1 >= r) return;
    int mid;
    mid = (l + r) / 2;
    mergesort(l , mid);
    mergesort(mid , r);

    merge(l , r , mid);
} 


-------------------------------------------------------------
                             ------- Graph ------
-------------------------------------------------------------
        bfs to find sortest move to x node to y node
-------------------------------------------------------------
bool mark[20][20];
int dis[30][30];
void bfs(int i, int j)
{
  for (int x = 1; x <= n; x++)
  {
    for (int y = 1; y <= m; y++)
    {
      dis[x][y] = INT_MAX;
      mark[x][y] = false;
    }
  }
  mark[i][j] = true;
  dis[i][j] = 0;
  queue<pair<int, int>>q;
  q.push(make_pair(i, j));
  while (!q.empty())
  {
    int x = q.front().first;
    int y = q.front().second;
    q.pop();
    for (int k = 0; k < 8; k++)
    {
      i = x + chessx[k];
      j = y + chessy[k];
      if (i >= 1 && i <= n && j >= 1 && j <= m && !mark[i][j])
      {
        dis[i][j] = min(dis[x][y] + 1, dis[i][j]);
        mark[i][j] = true;
        q.push(make_pair(i, j));
      }
    }
  }
}

-------------------------------------------------------------------
                            Data Stucture 
-------------------------------------------------------------------
      Segment Tree
-------------------------------------------------------------------
const int maxN = 1e5 + 7;
int ar[maxN];
int tree[4 * maxN];

void build(int at , int left , int right){
  if(left == right){
    tree[at] = ar[right];
    return ;
  }
  int mid = (left + right) >> 1 ;
  build(at * 2 , left , mid);
  build(at * 2 + 1 , mid + 1 , right);
  tree[at] = min(tree[at * 2] , tree[at * 2 + 1]);
}

int queries(int at , int left , int right , int l , int r){
  if(r < left or right < l)
    return (int)1e5 + 7 ;
  if(l <= left and right <= r)
    return tree[at];
  int mid = (left + right) >> 1 ;
  int x = queries(at * 2 , left , mid , l , r);
  int y = queries(at * 2 + 1 , mid + 1 , right , l , r);
  return min(x , y);
}

-------------------------------------------------------------------
                            DP
-------------------------------------------------------------------

-------------------------------------------------------------------
/*A Space Optimized DP solution for 0 - 1 Knapsack Problem*/
-------------------------------------------------------------------

ll knapSack(ll W, ll wt[], ll val[], ll n)
{
  ll dp[W + 1];
  memset(dp , 0 , sizeof(dp));

  for (int i = 0 ; i < n ; i ++) {
    for (int j = W ; j >= wt[i] ; j --) {
      dp[j] = max(dp[j] , val[i] + dp[j - wt[i]]);
    }
  }

  return dp[W];
}
-------------------------------------------------------------------
                         //sparse table
-------------------------------------------------------------------
const int MAXN = 1e5 + 1;

int logg[MAXN];
int ar[MAXN];
int st[MAXN][50];

void Logs() {
  logg[1] = 0 ;
  for (int i = 2 ; i < MAXN ; i ++)
    logg[i] = logg[i / 2] + 1;
}

void build(int n) {
  for (int i = 0 ; i < n ; i ++) {
    st[i][0] = ar[i];
  }

  for (int j = 1 ; j <= logg[n] ; j ++) {
    for (int i = 0 ; i + (1 << j) <= n ; i ++) {
      st[i][j] = min(st[i][j - 1] , st[i + (1 << (j - 1))][j - 1]);
    }
  }
}

int take_min(int a , int b) {
  int k = logg[b - a + 1];
  int ans = min(st[a][k] , st[b - (1 << k) + 1][k]);
  return ans ;
}
