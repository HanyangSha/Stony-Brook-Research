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
const ll inf = 1000000007;
 
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
int global_time;

int startTime[mxn];
int low[mxn];
int inverse_startTime[mxn];

vector<vi> adj;

vi st;
void build_st() {
    st.pb(root);
    F0R(i,n) if (i != root) {
        int l = st.size();
        int ridx = rand() % l;
        adj[st[ridx]].pb(i);
        st.pb(i);
    }
}

void tarjan(int v) {
    startTime[v] = global_time;
    low[v] = startTime[v];
    inverse_startTime[global_time] = v;
    global_time++;

    trav(x,adj[v]) {
        if (startTime[x] == -1) {
            tarjan(x);
            low[v] = min(low[v], low[x]);
        }
        low[v] = min(low[v], startTime[x]);
    }

    if (low[v] == startTime[v] && low[v] != 0) {
        int a = rand() % startTime[v]; // smaller start time: [0, startTime[v])
        int b = rand() % (global_time - startTime[v]) + startTime[v]; // same or greater start time: [startTime[v], global_time)
        int x = inverse_startTime[a];
        int y = inverse_startTime[b];
        adj[y].pb(x);
        low[v] = a;
    }
}

vector<pii> edges;
void add_edges(Eigen::MatrixXd &G) {
    F0R(i,n) F0R(j,n) if (i != j) if (G(i,j) == 0) { // edge does not exist 

        int p = rand() % 2; // 0 to 1
        if (p == 1) {
            adj[i].pb(j);
        }
        else {
            edges.pb(mp(i,j));
        }
    }
}

void init() {
    global_time = 0;
    memset(startTime, -1, sizeof startTime);
    memset(low, 0, sizeof low);
    memset(inverse_startTime, 0, sizeof inverse_startTime);
    adj.clear();
    st.clear();
    edges.clear();
}

double lambda2(Eigen::MatrixXd &L) {
    Eigen::EigenSolver<Eigen::MatrixXd> es;
    es.compute(L, /* computeEigenvectors = */ false);
    auto e = es.eigenvalues(); // returns a column vector of type complex<double>

    vector<double> eigenval;
    trav(x, e) eigenval.pb(x.real());
    sort(all(eigenval));

    return eigenval[1];
}

int main() {

    srand (time(NULL));

    bool check = true;

    FOR(n,5,6) { //F0R(iteration,50) {
        init();
        ::n = n;

        adj.resize(n);
        root = 0;
        build_st();

        tarjan(root); 

        Eigen::MatrixXd L(n,n);
        F0R(i,n) F0R(j,n) L(i,j) = 0;

        F0R(i,n) {
            trav(x,adj[i]) L(i,x) = 1;
        }

        add_edges(L);

        F0R(i,n) {
            //trav(x,adj[i]) L(i,x) = -1.0/adj[i].size();
            trav(x,adj[i]) L(i,x) = -1;
        }

        F0R(i,n) {
            cout << i << ": "; 
            trav(x,adj[i]) {
                cout << x << ' ';
            }
            cout << '\n';
        }

        F0R(i,n) L(i,i) = adj[i].size();

        cout << L << endl;

        double l1 = lambda2(L);

        // add 1 more edge

        if (edges.size() == 0) continue;

        int idx = rand() % edges.size();
        pii tmp = edges[idx];

        adj[tmp.first].pb(tmp.second);
        // L(tmp.first, tmp.second) = -1.0/adj[tmp.first].size();
        L(tmp.first, tmp.second) = -1;
        L(tmp.first, tmp.first)++;

        // trav(x, adj[tmp.first]) {
        //     L(tmp.first, x) = -1.0/adj[tmp.first].size();
        // }

        double l2 = lambda2(L);
        
        if (l2 - l1 < -1e-10) {
            check = false;
        }

        cout << L << endl;
        cout << check << endl;
    } //}

    return 0;
}