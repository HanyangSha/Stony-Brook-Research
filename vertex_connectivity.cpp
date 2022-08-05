#include <bits/stdc++.h>
#include "Eigen/Eigenvalues"
using namespace std;
 
typedef long long ll;
typedef vector<int> vi;
typedef vector<ll> vll;
typedef pair<int, int> pii;
typedef pair<ll, ll> pll;
 
#define FOR(i, a, b) for (int i=a; i<(b); i++)
#define F0R(i, a) for (int i=0; i<(a); i++)
#define FORd(i,a,b) for (int i = b; i >= a; i--)
#define F0Rd(i,a) for (int i = a; i >= 0; i--)
#define FORit(it,a) for (auto it = a.begin(); it != a.end(); it++)
#define trav(a,x) for (auto& a: x)
 
#define sz(x) (int)(x).size()
#define mp make_pair
#define pb push_back
#define all(x) x.begin(), x.end()
#define alla(arr, sz) arr, arr + sz
 
const int dx[4] = { 1,0,-1,0 }, dy[4] = { 0,1,0,-1 }; // for every grid problem!!
const ll linf = 4*1e18;
const int inf = 1000000007;
 
namespace io {
	void setIn(string s) { freopen(s.c_str(),"r",stdin); }
	void setOut(string s) { freopen(s.c_str(),"w",stdout); }
	void setIO(string s = "") {
		ios_base::sync_with_stdio(0); cin.tie(0); // fast I/O
		if (sz(s)) { setIn(s+".in"), setOut(s+".out"); } // for USACO
	}
}
 
using namespace io;

const int mxn = 1e5+5;

int n, root; 

vector<vi> adj;

vi st;
void build_st() {
    st.pb(root);
    F0R(i,n) if (i != root) {
        int l = st.size();
        int ridx = rand() % l;
        adj[st[ridx]].pb(i);
        adj[i].pb(st[ridx]);
        st.pb(i);
    }
}

void add_edges() {
    F0R(i,n) F0R(j,n) if (i != j) {
        bool b = true;
        trav(x, adj[i]) {
            if (x == j) {
                b = false; continue;
            }
        }
        if (b) {
            int p = rand() % 2; // 0 to 1
            if (p == 1) {
                adj[i].pb(j);
                adj[j].pb(i);
            }
        }
    }
}

double h2(int bd, int sd, int k) {
    double a = (k-1) * n * bd;
    double b = (n-k+1) * (k-1) + 4*(sd-k+2) * (n-sd-1);
    return a/b;
}

set<int> removed;
bool visited[20];

void dfs(int v) {
    visited[v] = true;
    trav(x, adj[v]) 
    if (removed.count(x) == 0) {
        if (!visited[x]) {
            dfs(x);
        }
    }
}

bool connect() {

    F0R(t,n) {
        if (removed.count(t) == 0 && !visited[t]) return false;
    }
    return true;
}

bool test(int k) {

    for (int i = 0; i <= (1<<n)-1; i++) 
    if (__builtin_popcount(i) == k) {

        removed.clear();
        memset(visited, false, sizeof visited);

        for(int j = 0; j <= n-1; j++) {
            if ((1<<j) & i) removed.insert(j);
        }

        F0R(j,n) {
            if (removed.count(j) == 0) {
                dfs(j);
                if (!connect()) return false;
                break;
            }
        }
    }

    return true; 
}

#define DEBUG 0

int main() {

    srand (time(NULL));

    #if !DEBUG

    adj.clear();
    st.clear();

    n = 5;

    adj.resize(n);
    root = 0;
    build_st();

    add_edges();

    #else

    setIn("test.in");
    int edge; 
    cin >> n >> edge;
    adj.resize(n);

    F0R(i,edge) {
        int a, b;
        cin >> a >> b;
        adj[a].pb(b);
        adj[b].pb(a);
    }

    #endif 

    Eigen::MatrixXd L(n,n);
    F0R(i,n) {
        F0R(j,n) L(i,j) = 0;
        L(i,i) = adj[i].size();
        trav(x,adj[i]) {
            L(i,x) = -1;
        }
    }

    Eigen::EigenSolver<Eigen::MatrixXd> es;
    es.compute(L, /* computeEigenvectors = */ false);
    auto e = es.eigenvalues(); // returns a column vector of type complex<double>

    vector<double> eigenval;
    trav(x, e) eigenval.pb(x.real());
    sort(all(eigenval));
    int l2 = eigenval[1];

    int sd = inf, bd = -1;
    F0R(i,n) {
        sd = min(sd, (int)(adj[i].size()));
        bd = max(bd, (int)(adj[i].size()));
    }

    int k = 2;
    for (; k <= sd; k++) {
        if (l2 <= h2(bd, sd, k)) break;
    }
    k--;
    
    cout << L << endl;
    cout << h2(bd, sd, k) << " " << l2 << " " << h2(bd, sd, k+1) << endl; // just use math, not related with actual connectivity after removing vertices
    cout << bd << " " << sd << " " << k << endl;

    cout << test(k-1) << endl; 
    cout << test(k);        
        

    return 0;
}