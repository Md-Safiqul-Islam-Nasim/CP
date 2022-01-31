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
--------------------------------------------------------------------------
bool cmp(pair < int , char > a , pair < int , char > b)
{
  if (a.first != b.first)
    return a.first > b.first;
  else
    return a.second < b.second ;
}
--------------------------------------------------------------------------
/*A Space Optimized DP solution for 0 - 1 Knapsack Problem*/

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
-------------------------------------------------------------------------- -
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
----------------------------------------------------------------------------
/*Decimal to Binary*/
vector < int > binaryNum ;
void decToBinary(int n)
{
  while (n > 0) {

    binaryNum.push_back(n % 2);
    n = n / 2;
  }
}

-------------------------------------------------------------------------- -
//convert hexadecimal to decimal
int convert(char num[]) {
  int len = strlen(num);
  int base = 1;
  int temp = 0;
  for (int i = len - 1; i >= 0; i--) {
    if (num[i] >= '0' && num[i] <= '9') {
      temp += (num[i] - 48) * base;
      base = base * 16;
    }
    else if (num[i] >= 'A' && num[i] <= 'F') {
      temp += (num[i] - 55) * base;
      base = base * 16;
    }
  }
  return temp;
}
-------------------------------------------------------------------------- -

int binaryToDecimal(int n)
{
  int num = n;
  int dec_value = 0;

  int base = 1;

  int temp = num;
  while (temp) {
    int last_digit = temp % 10;
    temp = temp / 10;

    dec_value += last_digit * base;

    base = base * 2;
  }

  return dec_value;
}
-------------------------------------------------------------------------- -
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

---------------------------------------------------------------------------- -
ll bigmod(ll a , ll p , ll m)
{
  if (p == 0)
    return 1; // If power is 0, then a ^ 0 = 1 for any value of a, And 1 Mod m=1 for any value of m, So return 1
  if (p % 2) // If power is odd, Split it : a ^ 5 😞 a )* (a ^ 4) --> left and right child respectively.
  {
    return ((a % m) * (bigmod(a , p - 1 , m))) % m;
  }
  else //If power is even then split it equally and return the result...
  {
    ll c = bigmod (a , p / 2 , m);
    return ((c % m) * (c % m)) % m;
  }
}
----------------------------------------------------------------------------
**If I have to face a large number of factorial then i have to use this formula**
firstly i have to find primefactor of that number

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
----------------------------------------------------------------------------
int binarySearch(int ar[], int x, int low, int high)
{

  while (low <= high)
  {
    int mid = low + (high - low) / 2;

    if (ar[mid] == x)
      return mid;

    if (ar[mid] < x)
      low = mid + 1;

    else
      high = mid - 1;
  }

  return -1;
}
------------------------------------------------------------------------------ -

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

------------------------------------------------------------------------------ -
int prime_factor (int num)
{
  int cnt = 0 ;
  int len = prime.size();

  for (int i = 0 ; i * i <= len ; i++)
  {
    if (num % prime[i] == 0)
    {
      while (num % prime[i] == 0)
      {
        cnt ++ ;

        num = num / prime[i];
      }

    }
  }

  if (num > 1)
  {
    cnt ++ ;
  }
  return cnt ;
}
-------------------------------------------------------------------------- -
int getSum(int n)
{
  int sum = 0;
  while (n != 0) {
    sum = sum + n % 10;
    n = n / 10;
  }
  return sum;
}
------------------------------------------------------------------------------ -
void simpleSieve(lld limit, vector<lld>& prime)
{
  bool mark[limit + 1];
  memset(mark, false, sizeof(mark));

  for (lld i = 2; i <= limit; ++i)
  {
    if (mark[i] == false)
    {
      prime.push_back(i);
      for (lld j = i; j <= limit; j += i)
        mark[j] = true;
    }
  }
}

void primesInRange(lld low, lld high)
{
  lld limit = floor(sqrt(high)) + 1;
  vector<lld> prime;
  simpleSieve(limit, prime);

  if (low == 1)
  {
    low++;
  }

  lld n = high - low + 1;

  bool mark[n + 1];
  memset(mark, false, sizeof(mark));

  for (lld i = 0; i < prime.size(); i++)
  {
    lld loLim = floor(low / prime[i]) * prime[i];
    if (loLim < low)
      loLim += prime[i];
    if (loLim == prime[i])
      loLim += prime[i];

    for (lld j = loLim; j <= high; j += prime[i])
      if (j != prime[i])
        mark[j - low] = true;
  }

  for (lld i = low; i <= high; i++)
    if (!mark[i - low])
      primeNumbers.push_back(i);
}

void callSieve(ull n)
{
  int range = n / 1000000;
  lld high = 1000000, low = 2;
  for (i = 0; i < range; i++)
  {
    primesInRange(low, high);
    low = high + 1;
    high += 1000000;
  }
}
---------------------------------------------------------------------------- -
**Factorial**
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
**Fibonacci**
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
------------------------------------------------------------ -
**Permutation Without Function**
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

--------------------------------------------------------------------
//NcR
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
------------------------------------------------------------------------
//Find NCR
ll p = r * (r - 1) ;
ll q = (r + b - 1) * (r + b);

*---------------------------------------------------------------------- -*
// Function to implement lower_bound
int lowerBound(int ara[], int low, int high, int key)
{
  int mid;
  while (low < high) {
    mid = (low + high) / 2;
    if (ara[mid] >= key)
      high = mid;
    else
      low = mid + 1;
  }
  return low;
}

int upperBound(int ara[], int low, int high, int key)
{
  int mid;
  while (low < high) {
    mid = (low + high) / 2;
    if (ara[mid] <= key)
      low = mid + 1;
    else
      high = mid;
  }
  return high;
}
---------------------------------------------------------------------------- -
//count co-prime of n
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
------------------------------------------------------------------------------
//count co-prime betwenn 1 to n
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
------------------------------------------------------------------------------
/*bfs to find sortest move to x node to y node*/
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